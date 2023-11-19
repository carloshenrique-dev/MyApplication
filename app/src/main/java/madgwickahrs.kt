import kotlin.math.sqrt

@Suppress("NAME_SHADOWING")
class MadgwickAHRS(private val sampleFreq: Float, private val gain: Float) {
    private var q0 = 1.0f
    private var q1: Float
    private var q2: Float
    private var q3 //quaternion
            = 0.0f

    init {
        q2 = q3
        q1 = q2
    }

    fun MadgwickAHRSupdate(
        gx: Float, gy: Float, gz: Float,
        ax: Float, ay: Float, az: Float,
        mx: Float, my: Float, mz: Float
    ) {
        var ax = ax
        var ay = ay
        var az = az
        var mx = mx
        var my = my
        var mz = mz
        var recipNorm: Float
        var s0: Float
        var s1: Float
        var s2: Float
        var s3: Float
        val hx: Float
        val hy: Float
        val _2q0mx: Float
        val _2q0my: Float
        val _2q0mz: Float
        val _2q1mx: Float
        val _2bx: Float
        val _2bz: Float
        val _4bx: Float
        val _4bz: Float
        val _2q0: Float
        val _2q1: Float
        val _2q2: Float
        val _2q3: Float
        val _2q0q2: Float
        val _2q2q3: Float
        val q0q0: Float
        val q0q1: Float
        val q0q2: Float
        val q0q3: Float
        val q1q1: Float
        val q1q2: Float
        val q1q3: Float
        val q2q2: Float
        val q2q3: Float
        val q3q3: Float

        // Use IMU algorithm if magnetometer measurement invalid (avoids NaN in magnetometer normalisation)
        if (mx == 0.0f && my == 0.0f && mz == 0.0f) {
            MadgwickAHRSupdateIMU(gx, gy, gz, ax, ay, az)
            return
        }

        // Rate of change of quaternion from gyroscope
        var qDot1: Float = 0.5f * (-q1 * gx - q2 * gy - q3 * gz)
        var qDot2: Float = 0.5f * (q0 * gx + q2 * gz - q3 * gy)
        var qDot3: Float = 0.5f * (q0 * gy - q1 * gz + q3 * gx)
        var qDot4: Float = 0.5f * (q0 * gz + q1 * gy - q2 * gx)

        // Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
        if (!(ax == 0.0f && ay == 0.0f && az == 0.0f)) {

            // Normalise accelerometer measurement
            recipNorm = invSqrt(ax * ax + ay * ay + az * az)
            ax *= recipNorm
            ay *= recipNorm
            az *= recipNorm

            // Normalise magnetometer measurement
            recipNorm = invSqrt(mx * mx + my * my + mz * mz)
            mx *= recipNorm
            my *= recipNorm
            mz *= recipNorm

            // Auxiliary variables to avoid repeated arithmetic
            _2q0mx = 2.0f * q0 * mx
            _2q0my = 2.0f * q0 * my
            _2q0mz = 2.0f * q0 * mz
            _2q1mx = 2.0f * q1 * mx
            _2q0 = 2.0f * q0
            _2q1 = 2.0f * q1
            _2q2 = 2.0f * q2
            _2q3 = 2.0f * q3
            _2q0q2 = 2.0f * q0 * q2
            _2q2q3 = 2.0f * q2 * q3
            q0q0 = q0 * q0
            q0q1 = q0 * q1
            q0q2 = q0 * q2
            q0q3 = q0 * q3
            q1q1 = q1 * q1
            q1q2 = q1 * q2
            q1q3 = q1 * q3
            q2q2 = q2 * q2
            q2q3 = q2 * q3
            q3q3 = q3 * q3

            // Reference direction of Earth's magnetic field
            hx =
                mx * q0q0 - _2q0my * q3 + _2q0mz * q2 + mx * q1q1 + _2q1 * my * q2 + _2q1 * mz * q3 - mx * q2q2 - mx * q3q3
            hy =
                _2q0mx * q3 + my * q0q0 - _2q0mz * q1 + _2q1mx * q2 - my * q1q1 + my * q2q2 + _2q2 * mz * q3 - my * q3q3
            _2bx = Math.sqrt((hx * hx + hy * hy).toDouble()).toFloat()
            _2bz =
                -_2q0mx * q2 + _2q0my * q1 + mz * q0q0 + _2q1mx * q3 - mz * q1q1 + _2q2 * my * q3 - mz * q2q2 + mz * q3q3
            _4bx = 2.0f * _2bx
            _4bz = 2.0f * _2bz

            // Gradient decent algorithm corrective step
            s0 =
                -_2q2 * (2.0f * q1q3 - _2q0q2 - ax) + _2q1 * (2.0f * q0q1 + _2q2q3 - ay) - _2bz * q2 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q3 + _2bz * q1) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q2 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz)
            s1 =
                _2q3 * (2.0f * q1q3 - _2q0q2 - ax) + _2q0 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * q1 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + _2bz * q3 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q2 + _2bz * q0) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q3 - _4bz * q1) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz)
            s2 =
                -_2q0 * (2.0f * q1q3 - _2q0q2 - ax) + _2q3 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * q2 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + (-_4bx * q2 - _2bz * q0) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q1 + _2bz * q3) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q0 - _4bz * q2) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz)
            s3 =
                _2q1 * (2.0f * q1q3 - _2q0q2 - ax) + _2q2 * (2.0f * q0q1 + _2q2q3 - ay) + (-_4bx * q3 + _2bz * q1) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q0 + _2bz * q2) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q1 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz)
            recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3) // normalise step magnitude
            s0 *= recipNorm
            s1 *= recipNorm
            s2 *= recipNorm
            s3 *= recipNorm

            // Apply feedback step
            qDot1 -= gain * s0
            qDot2 -= gain * s1
            qDot3 -= gain * s2
            qDot4 -= gain * s3
        }

        // Integrate rate of change of quaternion to yield quaternion
        q0 += qDot1 * (1.0f / sampleFreq)
        q1 += qDot2 * (1.0f / sampleFreq)
        q2 += qDot3 * (1.0f / sampleFreq)
        q3 += qDot4 * (1.0f / sampleFreq)

        // Normalise quaternion
        recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3)
        q0 *= recipNorm
        q1 *= recipNorm
        q2 *= recipNorm
        q3 *= recipNorm
    }

    fun MadgwickAHRSupdateIMU(gx: Float, gy: Float, gz: Float, ax: Float, ay: Float, az: Float) {
        var ax = ax
        var ay = ay
        var az = az
        var recipNorm: Float
        var s0: Float
        var s1: Float
        var s2: Float
        var s3: Float
        val _2q0: Float
        val _2q1: Float
        val _2q2: Float
        val _2q3: Float
        val _4q0: Float
        val _4q1: Float
        val _4q2: Float
        val _8q1: Float
        val _8q2: Float
        val q0q0: Float
        val q1q1: Float
        val q2q2: Float
        val q3q3: Float

        // Rate of change of quaternion from gyroscope
        var qDot1: Float = 0.5f * (-q1 * gx - q2 * gy - q3 * gz)
        var qDot2: Float = 0.5f * (q0 * gx + q2 * gz - q3 * gy)
        var qDot3: Float = 0.5f * (q0 * gy - q1 * gz + q3 * gx)
        var qDot4: Float = 0.5f * (q0 * gz + q1 * gy - q2 * gx)

        // Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
        if (!(ax == 0.0f && ay == 0.0f && az == 0.0f)) {

            // Normalise accelerometer measurement
            recipNorm = invSqrt(ax * ax + ay * ay + az * az)
            ax *= recipNorm
            ay *= recipNorm
            az *= recipNorm

            // Auxiliary variables to avoid repeated arithmetic
            _2q0 = 2.0f * q0
            _2q1 = 2.0f * q1
            _2q2 = 2.0f * q2
            _2q3 = 2.0f * q3
            _4q0 = 4.0f * q0
            _4q1 = 4.0f * q1
            _4q2 = 4.0f * q2
            _8q1 = 8.0f * q1
            _8q2 = 8.0f * q2
            q0q0 = q0 * q0
            q1q1 = q1 * q1
            q2q2 = q2 * q2
            q3q3 = q3 * q3

            // Gradient decent algorithm corrective step
            s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay
            s1 =
                _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az
            s2 =
                4.0f * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az
            s3 = 4.0f * q1q1 * q3 - _2q1 * ax + 4.0f * q2q2 * q3 - _2q2 * ay
            recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3) // normalise step magnitude
            s0 *= recipNorm
            s1 *= recipNorm
            s2 *= recipNorm
            s3 *= recipNorm

            // Apply feedback step
            qDot1 -= gain * s0
            qDot2 -= gain * s1
            qDot3 -= gain * s2
            qDot4 -= gain * s3
        }

        // Integrate rate of change of quaternion to yield quaternion
        q0 += qDot1 * (1.0f / sampleFreq)
        q1 += qDot2 * (1.0f / sampleFreq)
        q2 += qDot3 * (1.0f / sampleFreq)
        q3 += qDot4 * (1.0f / sampleFreq)

        // Normalise quaternion
        recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3)
        q0 *= recipNorm
        q1 *= recipNorm
        q2 *= recipNorm
        q3 *= recipNorm
    }

    companion object {
        private fun invSqrt(x: Float): Float {
            return (1.0f / sqrt(x.toDouble())).toFloat()
        }
    }
}