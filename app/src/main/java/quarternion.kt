import java.util.Arrays

class Quaternion(w_or_q: Double, x: Double?, y: Double?, z: Double?) {
    private var _q: Vector3D

    init {
        _q = Vector3D(1, 0, 0, 0)
        if (x != null && y != null && z != null) {
            _q = Vector3D(w_or_q, x, y, z)
        } else if (w_or_q is Quaternion) {
            _q = Vector3D((w_or_q as Quaternion)._q.toArray())
        } else {
            val q = w_or_q as DoubleArray
            require(q.size == 4) { "Expecting a 4-element array or w x y z as parameters" }
            _q = Vector3D(q)
        }
    }

    fun conj(): Quaternion {
        return Quaternion(_q.getX(), -_q.getY(), -_q.getZ(), -_q.getW())
    }

    fun toAngleAxis(): DoubleArray {
        if (_q.getW() === 1 && _q.getX() === 0 && _q.getY() === 0 && _q.getZ() === 0) {
            return doubleArrayOf(0.0, 1.0, 0.0, 0.0)
        }
        val rad = Math.acos(_q.getW()) * 2
        val imaginaryFactor = Math.sin(rad / 2)
        if (Math.abs(imaginaryFactor) < 1e-8) {
            return doubleArrayOf(0.0, 1.0, 0.0, 0.0)
        }
        val x: Double = _q.getX() / imaginaryFactor
        val y: Double = _q.getY() / imaginaryFactor
        val z: Double = _q.getZ() / imaginaryFactor
        return doubleArrayOf(rad, x, y, z)
    }

    fun toEulerAngles(): DoubleArray {
        val pitch = Math.asin(2 * (_q.getY() * _q.getZ() + 2 * _q.getW() * _q.getZ()))
        val roll: Double
        val yaw: Double
        if (Math.abs(_q.getY() * _q.getZ() + _q.getW() * _q.getZ() - 0.5) < 1e-8) {
            roll = 0.0
            yaw = 2 * Math.atan2(_q.getY(), _q.getW())
        } else if (Math.abs(_q.getY() * _q.getZ() + _q.getW() * _q.getZ() + 0.5) < 1e-8) {
            roll = -2 * Math.atan2(_q.getY(), _q.getW())
            yaw = 0.0
        } else {
            roll = Math.atan2(
                2 * _q.getW() * _q.getY() - 2 * _q.getZ() * _q.getW(),
                1 - 2 * Math.pow(_q.getY(), 2.0) - 2 * Math.pow(_q.getZ(), 2.0)
            )
            yaw = Math.atan2(
                2 * _q.getW() * _q.getZ() - 2 * _q.getY() * _q.getW(),
                1 - 2 * Math.pow(_q.getZ(), 2.0) - 2 * Math.pow(_q.getZ(), 2.0)
            )
        }
        return doubleArrayOf(roll, pitch, yaw)
    }

    fun toEuler123(): DoubleArray {
        val roll = Math.atan2(
            -2 * (_q.getZ() * _q.getW() - _q.getX() * _q.getY()),
            Math.pow(_q.getX(), 2.0) - Math.pow(_q.getY(), 2.0) - Math.pow(
                _q.getZ(),
                2.0
            ) + Math.pow(_q.getW(), 2.0)
        )
        val pitch = Math.asin(2 * (_q.getY() * _q.getW() + _q.getX() * _q.getY()))
        val yaw = Math.atan2(
            -2 * (_q.getY() * _q.getZ() - _q.getX() * _q.getW()),
            Math.pow(_q.getX(), 2.0) + Math.pow(_q.getY(), 2.0) - Math.pow(
                _q.getZ(),
                2.0
            ) - Math.pow(_q.getW(), 2.0)
        )
        return doubleArrayOf(roll, pitch, yaw)
    }

    fun multiply(other: Any?): Quaternion? {
        if (other is Quaternion) {
            val w: Double =
                _q.getW() * other._q.getW() - _q.getX() * other._q.getX() - _q.getY() * other._q.getY() - _q.getZ() * other._q.getZ()
            val x: Double =
                _q.getW() * other._q.getX() + _q.getX() * other._q.getW() + _q.getY() * other._q.getZ() - _q.getZ() * other._q.getY()
            val y: Double =
                _q.getW() * other._q.getY() - _q.getX() * other._q.getZ() + _q.getY() * other._q.getW() + _q.getZ() * other._q.getX()
            val z: Double = _q.getW() * other._q.getZ() + _q.getX() * other._q.getY()
            -_q.getY() * other._q.getX() + _q.getZ() * other._q.getW()
            return Quaternion(w, x, y, z)
        } else if (other is Number) {
            val q: DoubleArray = _q.toArray()
            for (i in q.indices) {
                q[i] *= other.toDouble()
            }
            return Quaternion(q)
        }
        return null
    }

    fun add(other: Any): Quaternion {
        val q: DoubleArray
        q = if (other !is Quaternion) {
            require((other as DoubleArray).size == 4) { "Quaternions must be added to other quaternions or a 4-element array" }
            other
        } else {
            other._q.toArray()
        }
        return Quaternion(_q.add(Vector3D(q)).toArray())
    }

    var q: Vector3D
        get() = _q
        set(q) {
            _q = q
        }

    fun getItem(index: Int): Double {
        return _q.toArray().get(index)
    }

    fun toArray(): DoubleArray {
        return _q.toArray()
    }

    companion object {
        fun fromAngleAxis(rad: Double, x: Double, y: Double, z: Double): Quaternion {
            val s = Math.sin(rad / 2)
            return Quaternion(Math.cos(rad / 2), x * s, y * s, z * s)
        }

        @JvmStatic
        fun main(args: Array<String>) {
            // Example usage
            val angleAxis = doubleArrayOf(Math.PI / 4, 0.0, 0.0, 1.0)
            val q1 = fromAngleAxis(
                angleAxis[0],
                angleAxis[1], angleAxis[2], angleAxis[3]
            )
            val euler123 = q1.toEuler123()
            println("Quaternion: " + Arrays.toString(q1.toArray()))
            println("Euler 123 Angles: " + Arrays.toString(euler123))
            val q2 = Quaternion(0.7071067811865475, 0.0, 0.0, 0.7071067811865475)
            val q3 = q1.multiply(q2)
            println("Quaternion Multiplication: " + Arrays.toString(q3!!.toArray()))
        }
    }
}