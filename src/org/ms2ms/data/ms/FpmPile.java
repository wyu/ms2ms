package org.ms2ms.data.ms;

public class FpmPile extends AbstractPile<FpmMatch>
{
  public FpmPile()       { super(); }
  public FpmPile(int s)  { super(s); init(); }

  @Override
  public int getKeyAt(int pile, int idx) { return get(pile,idx).getProtein(); }
  @Override
  public FpmMatch[] newPile(int s) { return new FpmMatch[s]; }

}