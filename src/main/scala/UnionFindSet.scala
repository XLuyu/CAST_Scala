import scala.collection.mutable.Map
/**
  * Created by workshop on 02-Nov-17.
  */
class UnionFindSet {
  private val father = Map[Int,Int]()
  def find(t: Int):Int = {
    if (!father.contains(t)) {father(t)=t; t} else iterFind(t)
  }
  def iterFind(t: Int):Int = {
    if (father(t)!=t) father(t) = iterFind(father(t))
    father(t)
  }

  def union(a: Int, b: Int) = {
    val A = find(a)
    val B = find(b)
    if (A!=B) father(Math.max(A,B)) = Math.min(A,B)
  }

  def isEmpty() = father.isEmpty
  def contains(elem:Int) = father.contains(elem)
  def getAnyKey() = father.last._1

  def printDebugInfo() = {
    for (entry <- father) {
      println(entry)
    }
  }
}
