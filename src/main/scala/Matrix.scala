/**
  * Created by workshop on 14-Apr-18.
  */
class Genotype(var acgt:Array[Double], val badCount:Int=0){
  val sum = acgt.sum
  if (sum!=0) acgt = acgt.map(_/sum)
  def isNonRepeat = (acgt.max>=0.9 || badCount<=0.1*sum) && badCount <= 0.9 * sum
  def isReliable = 5 <= sum && (acgt.max>=0.9 || badCount<=0.1*sum) && badCount <= 0.9 * sum //todo: 1<=sum
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
    var max = (0.0,0.0,0.0,0.0)
    var min = (1.0,1.0,1.0,1.0)
    for ( gt <- vector) {
      max = (Math.max(max._1,gt.acgt(0)),Math.max(max._2,gt.acgt(1)),Math.max(max._3,gt.acgt(2)),Math.max(max._4,gt.acgt(3)))
      min = (Math.min(min._1,gt.acgt(0)),Math.min(min._2,gt.acgt(1)),Math.min(min._3,gt.acgt(2)),Math.min(min._4,gt.acgt(3)))
    }
    (max._1-min._1)+(max._2-min._2)+(max._3-min._3)+(max._4-min._4)>=0.4
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
