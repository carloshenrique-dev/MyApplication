import kotlin.math.PI

fun normalizeAngles(angles: DoubleArray): DoubleArray {
    return angles.map { (it + PI) % (2 * PI) - PI }.toDoubleArray()
}