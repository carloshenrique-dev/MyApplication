import numpy as np
import sys
import json
from flask import Flask, request, jsonify

# Create a Flask app
app = Flask(__name__)

np.random.seed(777)

sys.path.append('./src')

from kalman_filters import ExtendedKalmanFilter as EKF
from utils import lla_to_enu, enu_to_lla ,normalize_angles, quaternion, madgwickahrs, savitzky_golay_filter

gt_trajectory_lla = []  # [longitude(deg), latitude(deg), altitude(meter)] x N
gt_yaws = []  # [yaw_angle(rad),] x N
gt_yaw_rates= []  # [vehicle_yaw_rate(rad/s),] x N
gt_forward_velocities = []  # [vehicle_forward_velocity(m/s),] x N
gt_timestamps = []
initialized = False
mad_filter = madgwickahrs.MadgwickAHRS(sampleperiod=1, quaternion=quaternion.Quaternion(1, 0, 0, 0), beta = 1, zeta = 0)

# Initial state estimate and covariance
x = np.array([0., 0., 0.])
P = np.array([
    [0., 0., 0.],
    [0., 0., 0.],
    [0., 0., 0.]
])

# Measurement error covariance Q
Q = np.array([
    [1., 0.],
    [0., 1.]
])

# State transition noise covariance R
r_value = 29.0
R = np.array([
    [r_value, 0., 0.],
    [0., r_value, 0.],
    [0., 0., r_value]
])

kf = EKF(x, P)

# Route to receive sensor data
@app.route('/sensor_data', methods=['POST'])
def receive_sensor_data():
    global  gt_timestamps, gt_trajectory_lla, gt_forward_velocities, gt_yaws, gt_yaw_rates ,kf, var_x, var_y, var_theta, mu_x, mu_y, mu_theta
    data = request.get_json()

    gt_trajectory_lla.append([
        data.get('latitude', 0), data.get('longitude', 0), data.get('altitude', 0)
    ])
    gt_forward_velocities.append(data.get('speed', 0))

    #Getting the Yaw with Madgwick AHRS filter
    
    gyroscope = [data.get('gyroscope', [0, 0, 0])]
    accelerometer = [data.get('accelerometer', [0, 0, 0])]
    #magnetometer = [df['x'][index], df['y'][index], df['z'][index]]
    mad_filter.update_imu(gyroscope, accelerometer)#, magnetometer)
    gt_timestamps = [data.get('timestamp', 0)]
    
    roll, pich, yaw = mad_filter.quaternion.to_euler_angles()

    if len(gt_yaws) > 0:
        last_yaw = gt_yaws[-1]
    else:
        # Set an initial value, e.g., 0 radians
        last_yaw = 0.0

    yaw_rate = yaw - last_yaw
    gt_yaws.append(yaw)
    gt_yaw_rates.append(yaw_rate)

    gt_trajectory_lla_new = np.array(gt_trajectory_lla).T
    gt_yaws_new = np.array(gt_yaws)
    gt_yaw_rates_new = np.array(gt_yaw_rates)
    gt_forward_velocities_new = np.array(gt_forward_velocities)

    lons, lats, _ = gt_trajectory_lla_new

    origin = gt_trajectory_lla_new[:, 0]  # set the initial position to the origin
    gt_trajectory_xyz = lla_to_enu(gt_trajectory_lla_new, origin)

    xs, ys, _ = gt_trajectory_xyz

    xy_obs_noise_std_1 = 2.0  # standard deviation of observation noise of x and y in meter

    xy_obs_noise = 1 #np.random.normal(0.0, xy_obs_noise_std_1, (2, ts))  # gen gaussian noise
    obs_trajectory_xyz1 = gt_trajectory_xyz.copy()
    obs_trajectory_xyz1[:2, :] += xy_obs_noise  # add the noise to ground-truth positions

    yaw_rate_noise_std = 0.2 # standard deviation of yaw rate in rad/s

    yaw_rate_noise = 1 #np.random.normal(0.0, yaw_rate_noise_std, (ts,))  # gen gaussian noise
    obs_yaw_rates = gt_yaw_rates_new.copy()
    obs_yaw_rates += yaw_rate_noise  # add the noise to ground-truth positions

    forward_velocity_noise_std = 1.0 # standard deviation of forward velocity in m/s

    forward_velocity_noise = 1 #np.random.normal(0.0, forward_velocity_noise_std, (ts,))  # gen gaussian noise
    obs_forward_velocities = gt_forward_velocities_new.copy()
    obs_forward_velocities += forward_velocity_noise  # add the noise to ground-truth positions

    initial_yaw_std = np.pi
    initial_yaw = gt_yaws_new[0] + np.random.normal(0, initial_yaw_std)
    xy_obs_noise_std = 0.119

    xy_obs_noise = 1 #np.random.normal(0.0, xy_obs_noise_std, (2, N)) 
    obs_trajectory_xyz = gt_trajectory_xyz.copy()
    obs_trajectory_xyz[:2, :] += xy_obs_noise  

    x = np.array([
        obs_trajectory_xyz[0, 0],
        obs_trajectory_xyz[1, 0],
        initial_yaw
    ])

    P = np.array([
        [xy_obs_noise_std ** 2., 0., 0.],
        [0., xy_obs_noise_std ** 2., 0.],
        [0., 0., initial_yaw_std ** 2.]
    ])


    Q = np.array([
        [xy_obs_noise_std ** 2., 0.],
        [0., xy_obs_noise_std ** 2.]
    ])


    R = np.array([
        [forward_velocity_noise_std ** 2., 0., 0.],
        [0., forward_velocity_noise_std ** 2., 0.],
        [0., 0., yaw_rate_noise_std ** 2.]
    ])


    global initialized
    if not initialized:
        # initialize Kalman filter
        kf = EKF(x, P)

        # array to store estimated 2d pose [x, y, theta]
        mu_x = [x[0],]
        mu_y = [x[1],]
        mu_theta = [x[2],]

        # array to store estimated error variance of 2d pose
        var_x = [P[0, 0],]
        var_y = [P[1, 1],]
        var_theta = [P[2, 2],]

    #NAO ESQUECER DE MUDAR ISSO
    dt = 1 #gt_timestamps[-1] - gt_timestamps[0]

    # get control input `u = [v, omega] + noise`
    u = np.array([
        obs_forward_velocities[-1],
        obs_yaw_rates[-1]
    ])
    
    # because velocity and yaw rate are multiplied with `dt` in state transition function,
    # its noise covariance must be multiplied with `dt**2.`
    R_ = R * (dt ** 2.)
    
    # propagate!
    kf.propagate(u, dt, R)
    
    # get measurement `z = [x, y] + noise`
    z = np.array([
        obs_trajectory_xyz[0, -1],
        obs_trajectory_xyz[1, -1]
    ])
    
    # update!
    kf.update(z, Q)
    
    # save estimated state to analyze later
    mu_x.append(kf.x[0])
    mu_y.append(kf.x[1])
    mu_theta.append(normalize_angles(kf.x[2]))
    
    # save estimated variance to analyze later
    var_x.append(kf.P[0, 0])
    var_y.append(kf.P[1, 1])
    var_theta.append(kf.P[2, 2])    

    mu_x_new = np.array(mu_x)
    mu_y_new = np.array(mu_y)
    mu_theta_new = np.array(mu_theta)

    var_x_new = np.array(var_x)
    var_y_new = np.array(var_y)
    var_theta_new = np.array(var_theta)

    xs, ys, _ = gt_trajectory_xyz
    xs, ys, _ = obs_trajectory_xyz1

    # Validate that mu_x and mu_y are of the same length
    if mu_x_new.shape[0] != mu_y_new.shape[0]:
        raise ValueError("mu_x and mu_y must be of the same length")

    # Create an array of zeros of the same shape as mu_x
    zeros_array = np.zeros(mu_x_new.shape)

    # Concatenate the arrays horizontally, this requires that the arrays be 2D,
    # so we use np.newaxis to increase their dimensions
    # Option 1: Using np.concatenate
    combined_array = np.concatenate((mu_x_new[:, np.newaxis], mu_y_new[:, np.newaxis], zeros_array[:, np.newaxis]), axis=1)

    combined_array.T

    lla_trajectory = enu_to_lla(combined_array.T,origin)

    lla_trajectory.T

    filtered_locations = [(loc[0],loc[1]) for loc in lla_trajectory.T]

    
    estimated_x, estimated_y, estimated_yaw = kf.x
    #result = enu_to_lla([estimated_x, estimated_y, 0], origin)
    #estimated_location = [float(x) for x in result.T[-1].tolist()]

    response = {
        'estimated_x': estimated_x,
        'estimated_y': estimated_y,
        'estimated_yaw': estimated_yaw,
        'estimated_location': filtered_locations[-1],
    }

    return jsonify(response)

if __name__ == "__main__":
    app.run('172.20.10.4')