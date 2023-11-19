import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.ArrayRealVector
import org.apache.commons.math3.linear.RealMatrix
import org.apache.commons.math3.linear.RealVector
import org.apache.commons.math3.linear.SingularMatrixException

class ExtendedKalmanFilter(x: RealVector, P: RealMatrix) {
    private var x // state to estimate: [x_, y_, theta]^T
            : RealVector
    private var P // estimation error covariance
            : RealMatrix

    init {
        this.x = x
        this.P = P
    }

    fun update(z: RealVector, Q: RealMatrix?) {
        // compute Kalman gain
        val H: RealMatrix = Array2DRowRealMatrix(
            arrayOf<DoubleArray>(
                doubleArrayOf(1.0, 0.0, 0.0),
                doubleArrayOf(0.0, 1.0, 0.0)
            )
        )
        val K: RealMatrix
        K = try {
            P.multiply(H.transpose())
                .multiply(H.multiply(P).multiply(H.transpose()).add(Q).inverse())
        } catch (e: SingularMatrixException) {
            throw RuntimeException("Singular matrix in update step")
        }

        // update state x
        val z_: RealVector = ArrayRealVector(doubleArrayOf(x.getEntry(0), x.getEntry(1)))
        x = x.add(K.operate(z.subtract(z_)))

        // update covariance P
        P = P.subtract(K.multiply(H).multiply(P))
    }

    fun propagate(u: RealVector, dt: Double, R: RealMatrix?) {
        // propagate state x
        val xVal: Double = x.getEntry(0)
        val yVal: Double = x.getEntry(1)
        val thetaVal: Double = x.getEntry(2)
        val v: Double = u.getEntry(0)
        val omega: Double = u.getEntry(1)
        val r = v / omega // turning radius
        val dtheta = omega * dt
        val dx = -r * Math.sin(thetaVal) + r * Math.sin(thetaVal + dtheta)
        val dy = r * Math.cos(thetaVal) - r * Math.cos(thetaVal + dtheta)
        x = x.add(ArrayRealVector(doubleArrayOf(dx, dy, dtheta)))

        // propagate covariance P
        val G: RealMatrix = Array2DRowRealMatrix(
            arrayOf<DoubleArray>(
                doubleArrayOf(
                    1.0,
                    0.0,
                    -r * Math.cos(thetaVal) + r * Math.cos(thetaVal + dtheta)
                ),
                doubleArrayOf(0.0, 1.0, -r * Math.sin(thetaVal) + r * Math.sin(thetaVal + dtheta)),
                doubleArrayOf(0.0, 0.0, 1.0)
            )
        ) // Jacobian of state transition function
        P = G.multiply(P).multiply(G.transpose()).add(R)
    }
}