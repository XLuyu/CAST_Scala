/**
  * Created by workshop on 14-Apr-18.
  */
class Genotype(var acgt:Array[Double], val badCount:Int=0){
  val sum = acgt.sum
  if (sum!=0) acgt = acgt.map(_/sum)
  def isNonRepeat = (acgt.max>=0.9 || badCount<=0.1*sum) && badCount <= 0.9 * sum
  def isReliable = 1 <= sum && (acgt.max>=0.9 || badCount<=0.1*sum) && badCount <= 0.9 * sum //todo: 1<=sum
  def distance(other:Genotype) = {
    if (this.isReliable && other.isReliable)
      (this.acgt zip other.acgt map {case (x,y)=>Math.abs(x-y)}).sum / 2
    else
      Double.NaN
  }
}

class GenotypeVector(val vector:Array[Genotype]){
  def toDistMatrix = {
      for (i <- vector) yield
        for (j <- vector) yield {
          val d = i.distance(j)
          if (d < 0.1) 0 else d
        }
  }
  def isReliable = vector.forall(_.isReliable)
  def isHeterogeneousPrecheck = {
    var max = Array(0.0,0.0,0.0,0.0,0.0)
    var min = Array(1.0,1.0,1.0,1.0,1.0)
    for ( gt <- vector) {
      max = max zip gt.acgt map {x=>Math.max(x._1,x._2)}
      min = min zip gt.acgt map {x=>Math.min(x._1,x._2)}
    }
    (max zip min map {x=>x._1-x._2}).sum>=0.4
  }
  def isHeterogeneous:Boolean = {
    if (!isHeterogeneousPrecheck) return false
    !toDistMatrix.exists(x=>{ val sx = x.sorted; sx.last-sx(1)<0.2}) // to check if heterogeneity is not from mutation
  }
  def consistentWithDepth(cov:Array[Double]) = {
    vector.indices.forall(i=>vector(i).sum<=cov(i)*1.8)

    /*todo:
      &&
      Util.PearsonCorrelationSimilarity(,cov)<0.3
      */
  }
}

class Matrix(var data:Array[Array[Double]]=null, size:Int=0){
  if (data==null) data = Array.ofDim[Double](size,size)
  def *=(factor:Double) = {
    for (i <- data.indices; j <- data(i).indices)
      data(i)(j) *= factor
    this
  }
  def +=(other:Matrix) = {
    for (i <- data.indices; j <- data(i).indices)
      data(i)(j) += other.data(i)(j)
    this
  }
  def copy = {
    val other = new Matrix(size=data.length)
    for ( i <- data.indices ; j <- data(i).indices)
      other.data(i)(j) = data(i)(j)
    other
  }
  def normalized = {
    val max = data.map(_.max)
    if (!max.contains(0.0)) {
      for ( i <- data.indices ; j <- data(i).indices)
        data(i)(j) /= max(i)
    }
    this
  }
  def printline(hint:String): Unit ={
    println(s"======$hint")
    for ( i <- data){
      i.foreach(j=>print(f"$j%.2f\t"))
      println()
    }
  }
}
