package org.ms2ms.alg;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;

/**
 * Created with IntelliJ IDEA.
 * User: hliu
 * Date: 8/19/14
 * Time: 5:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class XMLs
{
  public static XMLStreamReader newReader(String fileName) throws FileNotFoundException, XMLStreamException
  {
    return XMLInputFactory.newInstance().createXMLStreamReader(fileName, new FileInputStream(fileName));
  }
  public static boolean isA(XMLStreamReader parser, int a, String name)
  {
    return parser.getAttributeLocalName(a).equalsIgnoreCase(name);
  }
  public static Integer getInt(XMLStreamReader parser, int a)
  {
    return Integer.parseInt(parser.getAttributeValue(a));
  }
}

