/**
  * Created by workshop on 12-Jun-18.
  */
object Util {
  def PearsonCorrelationSimilarity(A: IndexedSeq[Double], B: IndexedSeq[Double]): Double = {
    val Asum = A.sum
    val Bsum = B.sum
    val A2sum = A.map(x => x * x).sum
    val B2sum = B.map(x => x * x).sum
    val product = A.zip(B).map(x => x._1 * x._2).sum
    val numerator = product - Asum * Bsum / A.length
    val dominator = math.sqrt((A2sum - Asum * Asum / A.length) * (B2sum - Bsum * Bsum / B.length))
    if (dominator == 0) 0 else numerator / dominator
  }
}
