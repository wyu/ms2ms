package org.ms2ms.utils;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class Timer
{
  //public static Debug debug = Debug.getInstance();
  private Map<String, StopWatch> watches = new ConcurrentHashMap<String, StopWatch>();

  public Timer() { super(); }

  public String start(String name) {
    // tack on the thread name for reference
    watches.put(threadedName(name), new StopWatch());
    return threadedName(name) + " started.";
  }
  public String stop(String name) {
    return threadedName(name) != null && watches.get(threadedName(name)) != null ? watches.get(threadedName(name)).stop() : "";
  }
  public Long msec(String name) {
    //return threadedName(name) + " stopped: " + watches.cells(threadedName(name)).stop();
    return threadedName(name) != null && watches.get(threadedName(name)) != null ? watches.get(threadedName(name)).msec() : 0L;
  }
  public String stop(String name, long base) {
    return threadedName(name) + " stopped: " + watches.get(threadedName(name)).stop(base);
  }
  public void pause ( String name) { watches.get(threadedName(name)).pause(); }
  public void resume (String name) { watches.get(threadedName(name)).resume(); }
  public String total(String name) {
    return threadedName(name) + " totaled: " + watches.get(threadedName(name)).total();
  }

  private String threadedName(String name) { return name + "#" + Thread.currentThread().getId(); }

  class StopWatch {
    protected Long mStart = 0L;
    protected Long mTotal = 0L;

    public StopWatch() { super(); start(); }

    public void   start()  { mStart  = System.currentTimeMillis(); mTotal = 0L; }
    public Long   msec()   { return    System.currentTimeMillis() - mStart; }
    public void   resume() { mStart  = System.currentTimeMillis(); }
    public void   pause()  { mTotal += System.currentTimeMillis() - mStart; }
    public String stop()   { return    Strs.msecToString(System.currentTimeMillis() - mStart); }
    public String total()  { return    Strs.msecToString(mTotal); }
    public String stop(long base)
    {
      return Strs.msecToString((System.currentTimeMillis() - mStart) / base);
    }
  }
}
