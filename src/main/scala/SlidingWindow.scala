
class SlidingWindowArray(val size:Int){
  val window = Array.fill(size)(Array(0,0,0,0,0))
  def inc(i:Int,j:Int) = window(i%size)(j) += 1
  def getAndClean(i:Int) ={
    val s = window(i%size)
    window(i%size) = Array(0,0,0,0,0)
    s
  }
}
class SlidingWindowInt(val size:Int){
  val window = Array.fill(size)(0)
  def inc(i:Int) = window(i%size) += 1
  def getAndClean(i:Int) ={
    val s = window(i%size)
    window(i%size) = 0
    s
  }
}