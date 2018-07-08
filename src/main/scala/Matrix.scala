/**
  * Created by workshop on 14-Apr-18.
  */
class Genotype(var acgt:Array[Double], val badCount:Int=0){
  val sum = acgt.sum
  if (sum!=0) acgt = acgt.map(_/sum)
//  def isNonRepeat = (acgt.max>=0.9 || badCount<=0.1*sum) && badCount <= 0.9 * sum
  val isReliable = badCount<=0.1*sum
  def distance(other:Genotype) = {
    if (this.isReliable && other.isReliable)
      if (this.sum==0 || other.sum==0) Double.NaN else (this.acgt zip other.acgt map {case (x,y)=>Math.abs(x-y)}).sum / 2
    else
      throw new Exception("[Error] try to compute unreliably genotypes' distance!!")
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
  val isReliable = vector.forall(_.isReliable) && vector.count(_.sum==0)<=vector.length/2
  def isHeterogeneousPrecheck = {
    var max = Array(0.0,0.0,0.0,0.0,0.0)
    var min = Array(1.0,1.0,1.0,1.0,1.0)
    for ( gt <- vector if gt.sum>0) {
      max = max zip gt.acgt map {x=>Math.max(x._1,x._2)}
      min = min zip gt.acgt map {x=>Math.min(x._1,x._2)}
    }
    (max zip min map {x=>x._1-x._2}).sum>=0.4
  }
  def isHeterogeneous:Boolean = {
    if (!isHeterogeneousPrecheck) return false
    !toDistMatrix.exists(x=>{ val sx = x.filter(!_.isNaN).sorted; sx.nonEmpty && sx.last-sx(1)<0.2}) // to check if heterogeneity is not from mutation
  }
  def consistentWithDepth(depth:Array[Double]) = {
    val fold = vector.indices.map(i=>vector(i).sum/depth(i))
    fold.forall(_<=1.8) &&
    toDistMatrix.map(dv=>Math.abs(Util.PearsonCorrelationSimilarity(dv,fold))).sum/vector.length<0.8
  }
}

class Matrix(var data:Array[Array[Double]]=null, size:Int=0){
  if (data==null) data = Array.ofDim[Double](size,size)
  def *=(factor:Double) = {
    for (i <- data.indices; j <- data(i).indices)
      data(i)(j) *= factor
    this
  }
  def *(factor:Double) = {
    val rt = this.copy
    for (i <- rt.data.indices; j <- rt.data(i).indices)
      rt.data(i)(j) *= factor
    rt
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
  def updateBy(other:Matrix): Unit ={
    for (i <- data.indices; j <- data(i).indices if !other.data(i)(j).isNaN)
      data(i)(j) = other.data(i)(j)
  }
  def fill(v:Double): Unit ={
    for (i <- data.indices; j <- data(i).indices) data(i)(j) = v
  }
  def containsNaN = data.exists(_.exists(_.isNaN))
  def printline(hint:String): Unit ={
    println(s"======$hint")
    for ( i <- data){
      i.foreach(j=>print(f"$j%.2f\t"))
      println()
    }
  }
}
