import kotlin.math.pow
import kotlin.math.absoluteValue

fun savitzkyGolay(y: DoubleArray, windowSize: Int, order: Int, deriv: Int = 0, rate: Int = 1): DoubleArray {
    val windowSizeAbs = windowSize.absoluteValue
    val orderAbs = order.absoluteValue

    if (windowSizeAbs % 2 != 1 || windowSizeAbs < 1) {
        throw IllegalArgumentException("windowSize must be a positive odd number")
    }

    if (windowSizeAbs < orderAbs + 2) {
        throw IllegalArgumentException("windowSize is too small for the polynomial order")
    }

    val orderRange = 0..orderAbs
    val halfWindow = (windowSizeAbs - 1) / 2

    // precompute coefficients
    val b = Array(windowSizeAbs) { k ->
        DoubleArray(orderAbs + 1) { i ->
            k.toDouble().pow(i.toDouble())
        }
    }

    val m = (b.pseudoInverse()[deriv].map { it * rate.toDouble().pow(deriv.toDouble()) * factorial(deriv) }).toDoubleArray()


    // pad the signal at the extremes with values taken from the signal itself
    val firstVals = y[0] - (y.sliceArray(1..halfWindow).reversedArray().map { it - y[0] }.maxOrNull()?.let { Math.abs(it.toDouble()) } ?: 0.0)
    val lastVals = y.last() + (y.sliceArray(y.size - halfWindow - 1 until y.size - 1).reversedArray().map { it - y.last() }.maxOrNull()?.let { Math.abs(it.toDouble()) } ?: 0.0)


    val extendedY = listOf(firstVals) + y.toList() + listOf(lastVals)

    return convolve(m.reversedArray(), extendedY)
}
fun factorial(n: Int): Double {
    var result = 1.0
    for (i in 2..n) {
        result *= i
    }
    return result
}
fun Array<DoubleArray>.pseudoInverse(): Array<DoubleArray> {
    val transpose = Array(size) { i -> DoubleArray(this[i].size) { j -> this[j][i] } }
    val product = Array(size) { i -> DoubleArray(this[i].size) { j -> this[i].zip(transpose[j]) { a, b -> a * b }.sum() } }
    return product.map { row -> row.map { it / size }.toDoubleArray() }.toTypedArray()
}

fun convolve(kernel: DoubleArray, signal: List<Double>): DoubleArray {
    val result = DoubleArray(signal.size - kernel.size + 1)
    for (i in result.indices) {
        result[i] = signal.subList(i, i + kernel.size).zip(kernel.toList()).map { it.first * it.second }.sum()
    }
    return result
}
