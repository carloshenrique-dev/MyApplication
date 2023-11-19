import java.util.Arrays
import kotlin.math.atan
import kotlin.math.cos
import kotlin.math.pow
import kotlin.math.sin
import kotlin.math.sqrt

object CoordinateTransformations {
    // Constant parameters defined in [1]
    private const val _a = 6378137.0
    private const val _f = 1.0 / 298.257223563
    private const val _b = (1.0 - _f) * _a
    private val _e = sqrt(_a * _a - _b * _b) / _a
    private val _e_prime = sqrt(_a * _a - _b * _b) / _b
    fun Rx(theta: Double): Array<DoubleArray> {
        val c = cos(theta)
        val s = sin(theta)
        return arrayOf(
            doubleArrayOf(1.0, 0.0, 0.0),
            doubleArrayOf(0.0, c, s),
            doubleArrayOf(0.0, -s, c)
        )
    }

    fun Ry(theta: Double): Array<DoubleArray> {
        val c = cos(theta)
        val s = sin(theta)
        return arrayOf(
            doubleArrayOf(c, 0.0, -s),
            doubleArrayOf(0.0, 1.0, 0.0),
            doubleArrayOf(s, 0.0, c)
        )
    }

    fun Rz(theta: Double): Array<DoubleArray> {
        val c = cos(theta)
        val s = sin(theta)
        return arrayOf(
            doubleArrayOf(c, s, 0.0),
            doubleArrayOf(-s, c, 0.0),
            doubleArrayOf(0.0, 0.0, 1.0)
        )
    }

    fun llaToEcef(pointsLla: DoubleArray): Array<DoubleArray> {
        val lon = Math.toRadians(pointsLla[0])
        val lat = Math.toRadians(pointsLla[1])
        val alt = pointsLla[2]
        val N = _a / sqrt(1.0 - (_e * sin(lat)).pow(2.0))
        val x = (N + alt) * cos(lat) * cos(lon)
        val y = (N + alt) * cos(lat) * sin(lon)
        val z = (N * (1.0 - _e.pow(2.0)) + alt) * sin(lat)
        return arrayOf(doubleArrayOf(x, y, z))
    }

    fun ecefToEnu(pointsEcef: Array<DoubleArray>, refLla: DoubleArray): Array<DoubleArray> {
        val lon = Math.toRadians(refLla[0])
        val lat = Math.toRadians(refLla[1])
        val refEcef = llaToEcef(refLla)
        val relative = doubleArrayOf(
            pointsEcef[0][0] - refEcef[0][0],
            pointsEcef[0][1] - refEcef[0][1],
            pointsEcef[0][2] - refEcef[0][2]
        )
        val R =
            multiplyMatrices(multiplyMatrices(Rz(Math.PI / 2.0), Ry(Math.PI / 2.0 - lat)), Rz(lon))
        val enuCoordinates = doubleArrayOf(
            R[0][0] * relative[0] + R[0][1] * relative[1] + R[0][2] * relative[2],
            R[1][0] * relative[0] + R[1][1] * relative[1] + R[1][2] * relative[2],
            R[2][0] * relative[0] + R[2][1] * relative[1] + R[2][2] * relative[2]
        )
        return arrayOf(
            doubleArrayOf(
                enuCoordinates[0],
                enuCoordinates[1], enuCoordinates[2]
            )
        )
    }

    fun llaToEnu(pointsLla: DoubleArray, refLla: DoubleArray): Array<DoubleArray> {
        val pointsEcef = llaToEcef(pointsLla)
        return ecefToEnu(pointsEcef, refLla)
    }

    fun enuToEcef(pointsEnu: Array<DoubleArray>, refLla: DoubleArray): Array<DoubleArray> {
        val lon = Math.toRadians(refLla[0])
        val lat = Math.toRadians(refLla[1])
        val refEcef = llaToEcef(refLla)
        val R =
            multiplyMatrices(multiplyMatrices(Rz(Math.PI / 2.0), Ry(Math.PI / 2.0 - lat)), Rz(lon))
        val relative = doubleArrayOf(
            R[0][0] * pointsEnu[0][0] + R[1][0] * pointsEnu[0][1] + R[2][0] * pointsEnu[0][2],
            R[0][1] * pointsEnu[0][0] + R[1][1] * pointsEnu[0][1] + R[2][1] * pointsEnu[0][2],
            R[0][2] * pointsEnu[0][0] + R[1][2] * pointsEnu[0][1] + R[2][2] * pointsEnu[0][2]
        )
        val ecefCoordinates = doubleArrayOf(
            refEcef[0][0] + relative[0],
            refEcef[0][1] + relative[1],
            refEcef[0][2] + relative[2]
        )
        return arrayOf(
            doubleArrayOf(
                ecefCoordinates[0],
                ecefCoordinates[1], ecefCoordinates[2]
            )
        )
    }

    fun ecefToLla(pointsEcef: Array<DoubleArray>): Array<DoubleArray> {
        val x = pointsEcef[0][0]
        val y = pointsEcef[0][1]
        val z = pointsEcef[0][2]
        val p = sqrt(x.pow(2.0) + y.pow(2.0))
        val theta = atan(z * _a / (p * _b))
        val lon = atan(y / x)
        val lat = atan(
            (z + _e_prime.pow(2.0) * _b * sin(theta).pow(3.0)) /
                    (p - _e.pow(2.0) * _a * cos(theta).pow(3.0))
        )
        val N = _a / sqrt(1.0 - (_e * sin(lat)).pow(2.0))
        val alt = p / cos(lat) - N
        return arrayOf(doubleArrayOf(Math.toDegrees(lon), Math.toDegrees(lat), alt))
    }

    fun enuToLla(pointsEnu: Array<DoubleArray>, refLla: DoubleArray): Array<DoubleArray> {
        val pointsEcef = enuToEcef(pointsEnu, refLla)
        return ecefToLla(pointsEcef)
    }

    private fun multiplyMatrices(a: Array<DoubleArray>, b: Array<DoubleArray>): Array<DoubleArray> {
        val aRows = a.size
        val aCols = a[0].size
        val bCols = b[0].size
        val result = Array(aRows) {
            DoubleArray(
                bCols
            )
        }
        for (i in 0 until aRows) {
            for (j in 0 until bCols) {
                for (k in 0 until aCols) {
                    result[i][j] += a[i][k] * b[k][j]
                }
            }
        }
        return result
    }

    @JvmStatic
    fun main(args: Array<String>) {
        // Example usage
        val refLla = doubleArrayOf(0.0, 0.0, 0.0)
        val pointLla = doubleArrayOf(1.0, 1.0, 0.0)
        val pointEnu = llaToEnu(pointLla, refLla)
        val pointLlaTransformed = enuToLla(pointEnu, refLla)
        println("Original Point (LLA): " + arrayOf(pointLla).contentDeepToString())
        println("Transformed Point (LLA): " + pointLlaTransformed.contentDeepToString())
    }
}