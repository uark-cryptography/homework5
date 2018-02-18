import java.math.BigInteger;
import java.util.Map;
import java.util.AbstractMap.SimpleImmutableEntry;
import java.util.ArrayList;
import java.util.Hashtable;

public class Cryptography {

    static private final String tab = "    ";
    static private final BigInteger ZERO = new BigInteger("0");

    // EXTENDED EUCLIDEAN

    public static BigInteger[] extendedEuclideanAlgorithm(BigInteger a, BigInteger b) {
        if (b.equals(BigInteger.ZERO)) {
            return new BigInteger[]{a, BigInteger.ONE, BigInteger.ZERO};
        }

        BigInteger u = BigInteger.ONE;
        BigInteger g = a;
        BigInteger x = BigInteger.ZERO;
        BigInteger y = b;

        while (!y.equals(BigInteger.ZERO)) {
            BigInteger t = g.mod(y);
            BigInteger q = g.subtract(t).divide(y);
            BigInteger s = u.subtract(q.multiply(x));

            u = x;
            g = y;
            x = s;
            y = t;
        }

        BigInteger v = (g.subtract(a.multiply(u))).divide(b);
        while (u.compareTo(BigInteger.ZERO) == -1) {
            u = u.add(b.divide(g));
            v = v.subtract(a.divide(g));
        }

        return new BigInteger[]{g, u, v};
    }

    static String formatExtendedEuclideanResults(BigInteger[] result, BigInteger a, BigInteger b) {
        return (
                result[1] + "·" + a + " + " +
                        result[2] + "·" + b + " = " +
                        result[0]
        );
    }

    static void runExtendedEuclidean() {
        BigInteger a = new BigInteger("2024");
        BigInteger b = new BigInteger("748");

        System.out.println("Example from book:");
        System.out.println(formatExtendedEuclideanResults(
                extendedEuclideanAlgorithm(a, b), a, b)
        );

        System.out.println("\n1.12 Exercise 'c':");

        System.out.print("  i)   ");
        a = new BigInteger("527");
        b = new BigInteger("1258");
        System.out.println(formatExtendedEuclideanResults(
                extendedEuclideanAlgorithm(a, b), a, b)
        );

        System.out.print("  ii)  ");
        a = new BigInteger("228");
        b = new BigInteger("1056");
        System.out.println(formatExtendedEuclideanResults(
                extendedEuclideanAlgorithm(a, b), a, b)
        );

        System.out.print("  iii) ");
        a = new BigInteger("163961");
        b = new BigInteger("167181");
        System.out.println(formatExtendedEuclideanResults(
                extendedEuclideanAlgorithm(a, b), a, b)
        );

        System.out.print("  iv)  ");
        a = new BigInteger("3892394");
        b = new BigInteger("239847");
        System.out.println(formatExtendedEuclideanResults(
                extendedEuclideanAlgorithm(a, b), a, b)
        );
    }

    // SHANK'S ALGORITHM

    // https://stackoverflow.com/questions/4407839/how-can-i-find-the-square-root-of-a-java-biginteger
    static BigInteger bigIntSqRootCeil(BigInteger x)
            throws IllegalArgumentException {
        if (x.compareTo(BigInteger.ZERO) < 0) {
            throw new IllegalArgumentException("Negative argument.");
        }

        if (x == BigInteger.ZERO || x == BigInteger.ONE) {
            return x;
        }

        BigInteger two = BigInteger.valueOf(2L);
        BigInteger y;
        for (
                y = x.divide(two);
                y.compareTo(x.divide(y)) > 0;
                y = ((x.divide(y)).add(y)).divide(two)
                )
            ;

        if (x.compareTo(y.multiply(y)) == 0) {
            return y;
        } else {
            return y.add(BigInteger.ONE);
        }
    }

    static BigInteger shanksAlgorithm(BigInteger g, BigInteger h, BigInteger N) {
        BigInteger n = bigIntSqRootCeil(N);

        Hashtable<BigInteger, BigInteger> table =
                new Hashtable<BigInteger, BigInteger>();

        for (
                BigInteger i = BigInteger.ZERO, key = BigInteger.ONE;
                i.compareTo(n) == -1;
                i = i.add(BigInteger.ONE),
                        key = key.multiply(g).mod(N)
                ) {
            table.put(key, i);
        }

        BigInteger u = g.modPow(n, N).modInverse(N);
        BigInteger value;
        for (
                BigInteger i = BigInteger.ZERO, key = h;
                i.compareTo(n) == -1;
                i = i.add(BigInteger.ONE),
                        key = key.multiply(u).mod(N)
                ) {
            value = table.get(key);
            if (value != null) {
                return i.multiply(n).add(value);
            }
        }

        throw new Error("No solution");
    }

    static void runShanksAlgorithm() {
        BigInteger g = new BigInteger("9704");
        BigInteger h = new BigInteger("13896");
        BigInteger N = new BigInteger("17389");

        System.out.println("Example from book:");
        System.out.println("9704^x ≅ 1389 (mod 17389)");
        System.out.println("x = " + shanksAlgorithm(g, h, N));

        String indent = "  ";
        System.out.println("\n2.17:");

        g = new BigInteger("11");
        h = new BigInteger("21");
        N = new BigInteger("71");
        System.out.println(indent + "a)");
        System.out.println(indent + indent + "11^x ≅ 21 (mod 71)");
        System.out.println(indent + indent + "x = " + shanksAlgorithm(g, h, N));

        g = new BigInteger("156");
        h = new BigInteger("116");
        N = new BigInteger("593");
        System.out.println(indent + "b)");
        System.out.println(indent + indent + "156^x ≅ 116 (mod 593)");
        System.out.println(indent + indent + "x = " + shanksAlgorithm(g, h, N));

        g = new BigInteger("650");
        h = new BigInteger("2213");
        N = new BigInteger("3571");
        System.out.println(indent + "b)");
        System.out.println(indent + indent + "650^x ≅ 2213 (mod 3571)");
        System.out.println(indent + indent + "x = " + shanksAlgorithm(g, h, N));

        g = new BigInteger("6");
        h = new BigInteger("7531");
        N = new BigInteger("8101");
        System.out.println("Other");
        System.out.println(indent + "6^x ≅ 7501 (mod 8101)");
        System.out.println(indent + "x = " + shanksAlgorithm(g, h, N));
    }

    // CHINESE REMAINDER THEROEM

    private static boolean areCoprime(ArrayList<BigInteger> list) {
        BigInteger iValue;
        BigInteger jValue;

        for(int i = 0; i < list.size(); i += 1) {
            iValue = list.get(i);
            for(int j = 0; j < list.size(); j += 1) {
                jValue = list.get(j);
                if (
                    !iValue.equals(jValue) &&
                    !iValue.gcd(jValue).equals(BigInteger.ONE)
                ) {
                    return false;
                }
            }
        }

        return true;
    }

    public static BigInteger chineseRemainderTheroem(
        ArrayList<SimpleImmutableEntry<BigInteger, BigInteger>> entries
    ) {
        BigInteger N = BigInteger.ONE;
        ArrayList<BigInteger> mods = new ArrayList<BigInteger>();
        for (int i = 0; i < entries.size(); i += 1) {
            BigInteger mod = entries.get(i).getValue(); 
            N = N.multiply(mod);
            mods.add(mod);
        };

        if (!areCoprime(mods)) {
            throw new Error("No solution - modulii are not coprime.");
        }

        ArrayList<BigInteger> inverse = new ArrayList<BigInteger>();
        BigInteger sum = ZERO;
        for (int i = 0; i < entries.size(); i += 1) {
            SimpleImmutableEntry<BigInteger, BigInteger> entry = entries.get(i);
            BigInteger value = entry.getKey();
            BigInteger mod = entry.getValue();

            BigInteger nextInverse = N.divide(mod).modInverse(mod);
            inverse.add(nextInverse);

            sum = sum.add(
                (N.divide(mod))
                    .multiply(value)
                    .multiply(nextInverse)
            );
        };

        if (sum.compareTo(ZERO) == -1) {
            return (sum.mod(N)).add(N);
        }
        return sum.mod(N);
    }

    public static void runChineseRemainderTheroem() {
        ArrayList<SimpleImmutableEntry<BigInteger, BigInteger>> entries = 
           new ArrayList<SimpleImmutableEntry<BigInteger, BigInteger>>();
        entries.add(new SimpleImmutableEntry<BigInteger, BigInteger>(
            new BigInteger("133"),
            new BigInteger("451")
        ));
        entries.add(new SimpleImmutableEntry<BigInteger, BigInteger>(
            new BigInteger("237"),
            new BigInteger("697")
        ));

        System.out.println("2.18");
        System.out.println(tab + "c)");
        try {
            chineseRemainderTheroem(entries);
        } catch (Error error) {
            System.out.println(tab + tab + error.getMessage());
        }

        entries.clear();
        entries.add(new SimpleImmutableEntry<BigInteger, BigInteger>(
            new BigInteger("37"),
            new BigInteger("43")
        ));
        entries.add(new SimpleImmutableEntry<BigInteger, BigInteger>(
            new BigInteger("22"),
            new BigInteger("49")
        ));
        entries.add(new SimpleImmutableEntry<BigInteger, BigInteger>(
            new BigInteger("18"),
            new BigInteger("71")
        ));
        System.out.println(tab + "e)");
        System.out.println(tab +  tab + chineseRemainderTheroem(entries));
    }

    // POHLIG-HELLMAN ALGORITHM

    private static BigInteger pow(BigInteger base, BigInteger exponent) {
        BigInteger result = BigInteger.ONE;
        while (exponent.signum() > 0) {
            if (exponent.testBit(0)) result = result.multiply(base);
            base = base.multiply(base);
            exponent = exponent.shiftRight(1);
        }
        return result;
    }

    private static ArrayList<BigInteger> primeFactors(BigInteger n) {
        BigInteger two = BigInteger.valueOf(2);
        ArrayList<BigInteger> fs = new ArrayList<BigInteger>();

        if (n.compareTo(two) < 0)
        {
            throw new IllegalArgumentException("must be greater than one");
        }

        while (n.mod(two).equals(BigInteger.ZERO))
        {
            fs.add(two);
            n = n.divide(two);
        }

        if (n.compareTo(BigInteger.ONE) > 0)
        {
            BigInteger f = BigInteger.valueOf(3);
            while (f.multiply(f).compareTo(n) <= 0)
            {
                if (n.mod(f).equals(BigInteger.ZERO))
                {
                    fs.add(f);
                    n = n.divide(f);
                }
                else
                {
                    f = f.add(two);
                }
            }
            fs.add(n);
        }

        return fs;
    }

    private static ArrayList<SimpleImmutableEntry<BigInteger, BigInteger>> toListOfTuples(
        ArrayList<BigInteger> list
    ) {
        Hashtable<BigInteger, BigInteger> table =
            new Hashtable<BigInteger, BigInteger>();

        list.forEach(key -> {
            BigInteger value = table.get(key);
            if (value == null) {
                table.put(key, BigInteger.ONE);
            } else {
                table.replace(key, value.add(BigInteger.ONE));
            }
        });

        ArrayList<SimpleImmutableEntry<BigInteger, BigInteger>> listOfTuple =
            new ArrayList<SimpleImmutableEntry<BigInteger, BigInteger>>();
        for(Map.Entry<BigInteger, BigInteger> entry : table.entrySet()) {
            BigInteger key = entry.getKey();
            BigInteger value = entry.getValue();
            listOfTuple.add(new SimpleImmutableEntry<BigInteger, BigInteger>(
                key,
                value
            ));
        }

        return listOfTuple;
    }

    public static BigInteger pohligHellmanAlgorithm(
        BigInteger g,
        BigInteger h,
        BigInteger N
    ) {

        ArrayList<SimpleImmutableEntry<BigInteger, BigInteger>> smallPrimes =
            toListOfTuples(primeFactors(N.subtract(BigInteger.ONE)));

        ArrayList<SimpleImmutableEntry<BigInteger, BigInteger>> values =
            new ArrayList<SimpleImmutableEntry<BigInteger, BigInteger>>();

        smallPrimes.forEach((tuple) -> {
            BigInteger base = tuple.getKey();
            BigInteger exp = tuple.getValue();
            BigInteger power = pow(base, exp);

            BigInteger p = pow(base, exp.subtract(BigInteger.ONE));
            BigInteger subG = pow(g, p).mod(N);
            BigInteger subH = pow(h, p).mod(N);
            BigInteger c;
            BigInteger cHelper = null;
            BigInteger x = null;
            for (
                BigInteger j = BigInteger.ONE,
                curPow = BigInteger.ONE;
                j.compareTo(exp) < 1;
                j = j.add(BigInteger.ONE),
                curPow = curPow.multiply(base)
            ) {
                c = shanksAlgorithm(subG, subH, N).mod(base);

                if (x == null) {
                    x = c;
                } else {
                    x = x.add(c.multiply(curPow));
                }

                p = pow(base, exp.subtract(j.add(BigInteger.ONE)));
                if (cHelper == null) {
                    cHelper = c.negate().multiply(curPow);
                } else {
                    cHelper = cHelper.add(c.negate().multiply(curPow));
                }
                subH = pow(h.multiply(g.modPow(cHelper, N)), p).mod(N);
            }
            values.add(new SimpleImmutableEntry<BigInteger, BigInteger>(
                x,
                power 
            ));
        });

        return chineseRemainderTheroem(values);
    }

    public static void runPohligHellmanAlgorithm() {
        BigInteger g = new BigInteger("6");
        BigInteger h = new BigInteger("7531");
        BigInteger N = new BigInteger("8101");

        System.out.println("Example from class");
        System.out.println(tab + g + "^x = " + h + " mod " + N);
        System.out.println(tab + "x = " + pohligHellmanAlgorithm(g, h, N));

        System.out.println("2.28");

        g = new BigInteger("7");
        h = new BigInteger("166");
        N = new BigInteger("443");
        System.out.println(tab + "a)");
        System.out.println(tab + tab + g + "^x = " + h + " mod " + N);
        System.out.println(tab + tab + "x = " + pohligHellmanAlgorithm(g, h, N));

        g = new BigInteger("10");
        h = new BigInteger("243278");
        N = new BigInteger("746497");
        System.out.println(tab + "b)");
        System.out.println(tab + tab + g + "^x = " + h + " mod " + N);
        System.out.println(tab + tab + "x = " + pohligHellmanAlgorithm(g, h, N));

        g = new BigInteger("2");
        h = new BigInteger("39183497");
        N = new BigInteger("41022299");
        System.out.println(tab + "c)");
        System.out.println(tab + tab + g + "^x = " + h + " mod " + N);
        System.out.println(tab + tab + "x = " + pohligHellmanAlgorithm(g, h, N));

        g = new BigInteger("17");
        h = new BigInteger("192988");
        N = new BigInteger("1291799");
        System.out.println(tab + "d)");
        System.out.println(tab + tab + g + "^x = " + h + " mod " + N);
        System.out.println(tab + tab + "x = " + pohligHellmanAlgorithm(g, h, N));
    }

    // RUNNER

    public static void main(String[] args) {
        // runExtendedEuclidean();
        // runShanksAlgorithm();
        // runChineseRemainderTheroem();
        runPohligHellmanAlgorithm();
    }
}
