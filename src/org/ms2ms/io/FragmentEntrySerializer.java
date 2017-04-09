package org.ms2ms.io;

import org.mapdb.DataInput2;
import org.mapdb.DataOutput2;
import org.mapdb.Serializer;
import org.ms2ms.data.ms.FragmentEntry;
import org.ms2ms.utils.IOs;

import java.io.IOException;
import java.util.Comparator;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.function.ToLongFunction;

/**
 * Created by yuw on 3/31/17.
 */
public class FragmentEntrySerializer implements Serializer<FragmentEntry>
{
  /**
   * Serializes the content of the given value into the given
   * {@link DataOutput2}.
   *
   * @param out   DataOutput2 to save object into
   * @param value Object to serialize
   * @throws IOException in case of an I/O error
   */
  @Override
  public void serialize(DataOutput2 out, FragmentEntry value) throws IOException
  {
    if (value==null)
      System.out.print("");
    out.writeInt(  value.getLength());
    out.writeByte( value.getCharge());
    out.writeInt(  value.getPeptideKey());
    out.writeFloat(value.getMH());
    // can't write the actual object. to avoid recursive
    // write the peptide seq key and MH so we can hook up the right frag at later time
    out.writeFloat(value.getPrev()!=null?value.getPrev().getMH():0f);
  }

  /**
   * Deserializes and returns the content of the given {@link DataInput2}.
   *
   * @param input     DataInput2 to de-serialize data from
   * @param available how many bytes that are available in the DataInput2 for
   *                  reading, may be -1 (in streams) or 0 (null).
   * @return the de-serialized content of the given {@link DataInput2}
   * @throws IOException in case of an I/O error
   */
  @Override
  public FragmentEntry deserialize(DataInput2 input, int available) throws IOException
  {
    if (available>0)
    {
      FragmentEntry F = new FragmentEntry();
      F.setLen(input.readInt()).setCharge(input.readChar()).setPeptideKey(input.readInt()).setMH(input.readFloat());

      float p = input.readFloat();
      if (p!=0f) F.setPrev(new FragmentEntry(p, '0', 0, null, 0));
    }
    return null;
  }

  /**
   * Returns the fixed size of the serialized form in bytes or -1 if the size
   * is not fixed (e.g. for Strings).
   * <p>
   * Some optimizations can be applied to serializers with a fixed size.
   *
   * @return the fixed size of the serialized form in bytes or -1 if the size
   * is not fixed
   */
  @Override
  public int fixedSize() { return 4+1+4+4+4; }

  /**
   * Returns if this Serializer is trusted to always read the same number of
   * bytes as it writes for any given object being serialized/de-serialized.
   * <p>
   * MapDB has a relaxed record size boundary checking. It expects
   * deserializers to read exactly as many bytes as were written during
   * serialization. If a deserializer reads more bytes than it wrote, it might
   * start reading others record data in store.
   * <p>
   * Some serializers (Kryo) have problems with this. To prevent this, we can
   * not read data directly from a store, but we must copy them into separate
   * {@code byte[]} buffers. Thus, zero-copy optimizations are disabled by
   * default, but can be explicitly enabled here by letting this method return
   * {@code true}.
   * <p>
   * This flag indicates if this serializer was 'verified' to read as many
   * bytes as it writes. It should also be much better tested etc.
   *
   * @return if this Serializer is trusted to always read the same number of
   * bytes as it writes for any given object being serialized/de-serialized
   */
  @Override
  public boolean isTrusted() { return true; }

  @Override
  public int compare(FragmentEntry first, FragmentEntry second) { return first.compareTo(second); }

  /**
   * Returns a comparator that imposes the reverse ordering of this
   * comparator.
   *
   * @return a comparator that imposes the reverse ordering of this
   * comparator.
   * @since 1.8
   */
  @Override
  public Comparator<FragmentEntry> reversed()
  {
    return null;
  }

  /**
   * Returns a lexicographic-order comparator with another comparator.
   * If this {@code Comparator} considers two elements equal, i.e.
   * {@code compare(a, b) == 0}, {@code other} is used to determine the order.
   * <p>
   * <p>The returned comparator is serializable if the specified comparator
   * is also serializable.
   *
   * @param other the other comparator to be used when this comparator
   *              compares two objects that are equal.
   * @return a lexicographic-order comparator composed of this and then the
   * other comparator
   * @throws NullPointerException if the argument is null.
   * @apiNote For example, to sort a collection of {@code String} based on the length
   * and then case-insensitive natural ordering, the comparator can be
   * composed using following code,
   * <p>
   * <pre>{@code
   *     Comparator<String> cmp = Comparator.comparingInt(String::length)
   *             .thenComparing(String.CASE_INSENSITIVE_ORDER);
   * }</pre>
   * @since 1.8
   */
  @Override
  public Comparator<FragmentEntry> thenComparing(Comparator<? super FragmentEntry> other)
  {
    return null;
  }

  /**
   * Returns a lexicographic-order comparator with a function that
   * extracts a key to be compared with the given {@code Comparator}.
   *
   * @param keyExtractor  the function used to extract the sort key
   * @param keyComparator the {@code Comparator} used to compare the sort key
   * @return a lexicographic-order comparator composed of this comparator
   * and then comparing on the key extracted by the keyExtractor function
   * @throws NullPointerException if either argument is null.
   * @implSpec This default implementation behaves as if {@code
   * thenComparing(comparing(keyExtractor, cmp))}.
   * @see #comparing(Function, Comparator)
   * @see #thenComparing(Comparator)
   * @since 1.8
   */
  @Override
  public <U> Comparator<FragmentEntry> thenComparing(Function<? super FragmentEntry, ? extends U> keyExtractor, Comparator<? super U> keyComparator)
  {
    return null;
  }

  /**
   * Returns a lexicographic-order comparator with a function that
   * extracts a {@code Comparable} sort key.
   *
   * @param keyExtractor the function used to extract the {@link
   *                     Comparable} sort key
   * @return a lexicographic-order comparator composed of this and then the
   * {@link Comparable} sort key.
   * @throws NullPointerException if the argument is null.
   * @implSpec This default implementation behaves as if {@code
   * thenComparing(comparing(keyExtractor))}.
   * @see #comparing(Function)
   * @see #thenComparing(Comparator)
   * @since 1.8
   */
  @Override
  public <U extends Comparable<? super U>> Comparator<FragmentEntry> thenComparing(Function<? super FragmentEntry, ? extends U> keyExtractor)
  {
    return null;
  }

  /**
   * Returns a lexicographic-order comparator with a function that
   * extracts a {@code int} sort key.
   *
   * @param keyExtractor the function used to extract the integer sort key
   * @return a lexicographic-order comparator composed of this and then the
   * {@code int} sort key
   * @throws NullPointerException if the argument is null.
   * @implSpec This default implementation behaves as if {@code
   * thenComparing(comparingInt(keyExtractor))}.
   * @see #comparingInt(ToIntFunction)
   * @see #thenComparing(Comparator)
   * @since 1.8
   */
  @Override
  public Comparator<FragmentEntry> thenComparingInt(ToIntFunction<? super FragmentEntry> keyExtractor)
  {
    return null;
  }

  /**
   * Returns a lexicographic-order comparator with a function that
   * extracts a {@code long} sort key.
   *
   * @param keyExtractor the function used to extract the long sort key
   * @return a lexicographic-order comparator composed of this and then the
   * {@code long} sort key
   * @throws NullPointerException if the argument is null.
   * @implSpec This default implementation behaves as if {@code
   * thenComparing(comparingLong(keyExtractor))}.
   * @see #comparingLong(ToLongFunction)
   * @see #thenComparing(Comparator)
   * @since 1.8
   */
  @Override
  public Comparator<FragmentEntry> thenComparingLong(ToLongFunction<? super FragmentEntry> keyExtractor)
  {
    return null;
  }

  /**
   * Returns a lexicographic-order comparator with a function that
   * extracts a {@code double} sort key.
   *
   * @param keyExtractor the function used to extract the double sort key
   * @return a lexicographic-order comparator composed of this and then the
   * {@code double} sort key
   * @throws NullPointerException if the argument is null.
   * @implSpec This default implementation behaves as if {@code
   * thenComparing(comparingDouble(keyExtractor))}.
   * @see #comparingDouble(ToDoubleFunction)
   * @see #thenComparing(Comparator)
   * @since 1.8
   */
  @Override
  public Comparator<FragmentEntry> thenComparingDouble(ToDoubleFunction<? super FragmentEntry> keyExtractor)
  {
    return null;
  }

  /**
   * Returns if the first and second arguments are equal to each other.
   * Consequently, if both arguments are {@code null}, {@code true} is
   * returned and if exactly one argument is {@code null}, {@code false} is
   * returned.
   *
   * @param first  an object
   * @param second another object to be compared with the first object for
   *               equality
   * @return if the first and second arguments are equal to each other
   * @see Object#equals(Object)
   */
  @Override
  public boolean equals(FragmentEntry first, FragmentEntry second) { return first.compareTo(second)==0; }

  /**
   * Returns a hash code of a given non-null argument. The output of the
   * method is affected by the given seed, allowing protection against crafted
   * hash attacks and to provide a better distribution of hashes.
   *
   * @param o    an object
   * @param seed used to "scramble" the
   * @return a hash code of a non-null argument
   * @throws NullPointerException if the provided object is null
   * @see Object#hashCode
   */
  @Override
  public int hashCode(FragmentEntry o, int seed) { return o.hashCode(); }

  /**
   * TODO: Document this method
   *
   * @return
   */
  @Override
  public boolean needsAvailableSizeHint() { return false; }

  /**
   * Deserializes and returns the content of the given long.
   *
   * @param input     long to de-serialize data from
   * @param available how many bytes that are available in the long for
   *                  reading, or 0 (null).
   * @return the de-serialized content of the given long
   * @throws IOException in case of an I/O error
   */
  @Override
  public FragmentEntry deserializeFromLong(long input, int available) throws IOException
  {
    return null;
  }

  /**
   * Creates binary copy of given object. If the datatype is immutable the same instance might be returned
   *
   * @param value
   */
  @Override
  public FragmentEntry clone(FragmentEntry value) throws IOException { return value!=null?new FragmentEntry(value):null; }
}
