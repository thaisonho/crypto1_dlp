#pragma GCC optimize("O3")
#include <iostream>
#include <fstream>
#include <string>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <vector>
#include <sstream>
#include <random>
#include <ctime>
#include <cstdlib>
#include <iomanip>

using namespace std;

class BigInt {
private:
    vector<uint32_t> data;

    void normalize() {
        while (data.size() > 1 && data.back() == 0) {
            data.pop_back();
        }
        if (data.empty()) {
            data.push_back(0);
        }
    }

public:
    BigInt() {
        data.push_back(0);
    }

    BigInt(uint64_t val) {
        data.push_back((uint32_t)(val & 0xFFFFFFFF));
        if (val > 0xFFFFFFFF) {
            data.push_back((uint32_t)(val >> 32));
        }
        normalize();
    }

    BigInt(const string& hex) {
        if (hex.empty() || hex == "0") {
            data.push_back(0);
            return;
        }

        data.reserve((hex.length() + 7) / 8);
        data.push_back(0);

        for (int i = 0; i < (int)hex.length(); i++) {
            char c = hex[i];
            int digit;
            if (c >= '0' && c <= '9') digit = c - '0';
            else if (c >= 'A' && c <= 'F') digit = c - 'A' + 10;
            else if (c >= 'a' && c <= 'f') digit = c - 'a' + 10;
            else continue;

            int bitPos = i * 4;
            int wordPos = bitPos / 32;
            int bitInWord = bitPos % 32;

            if (wordPos >= (int)data.size()) {
                data.resize(wordPos + 1, 0);
            }

            data[wordPos] |= ((uint32_t)digit << bitInWord);
            if (bitInWord > 28) {
                if (wordPos + 1 >= (int)data.size()) {
                    data.resize(wordPos + 2, 0);
                }
                data[wordPos + 1] |= ((uint32_t)digit >> (32 - bitInWord));
            }
        }

        normalize();
    }

    bool isZero() const { return data.size() == 1 && data[0] == 0; }
    bool isOne() const { return data.size() == 1 && data[0] == 1; }
    bool isEven() const { return (data[0] & 1) == 0; }

    bool getBit(int pos) const {
        int wordPos = pos / 32;
        int bitPos = pos % 32;
        if (wordPos >= (int)data.size()) return false;
        return (data[wordPos] >> bitPos) & 1;
    }

    void setBit(int pos) {
        int wordPos = pos / 32;
        int bitPos = pos % 32;
        if (wordPos >= (int)data.size()) {
            data.resize(wordPos + 1, 0);
        }
        data[wordPos] |= (1U << bitPos);
    }

    int bitLength() const {
        if (isZero()) return 0;
        int len = (data.size() - 1) * 32;
        uint32_t top = data.back();
        while (top > 0) { len++; top >>= 1; }
        return len;
    }

    bool operator==(const BigInt& other) const { return data == other.data; }
    bool operator!=(const BigInt& other) const { return !(*this == other); }

    bool operator<(const BigInt& other) const {
        if (data.size() != other.data.size()) return data.size() < other.data.size();
        for (int i = data.size() - 1; i >= 0; i--) {
            if (data[i] != other.data[i]) return data[i] < other.data[i];
        }
        return false;
    }

    bool operator<=(const BigInt& other) const { return !(other < *this); }
    bool operator>(const BigInt& other) const { return other < *this; }
    bool operator>=(const BigInt& other) const { return !(*this < other); }

    BigInt operator+(const BigInt& other) const {
        BigInt result;
        result.data.clear();
        uint64_t carry = 0;
        int maxSize = max(data.size(), other.data.size());
        result.data.reserve(maxSize + 1);

        for (int i = 0; i < maxSize || carry; i++) {
            uint64_t sum = carry;
            if (i < (int)data.size()) sum += data[i];
            if (i < (int)other.data.size()) sum += other.data[i];
            result.data.push_back((uint32_t)(sum & 0xFFFFFFFF));
            carry = sum >> 32;
        }
        result.normalize();
        return result;
    }

    BigInt operator-(const BigInt& other) const {
        if (*this < other) return BigInt(0);
        BigInt result;
        result.data.resize(data.size());
        int64_t borrow = 0;

        for (int i = 0; i < (int)data.size(); i++) {
            int64_t diff = (int64_t)data[i] - borrow;
            if (i < (int)other.data.size()) diff -= other.data[i];
            if (diff < 0) { diff += 0x100000000LL; borrow = 1; }
            else { borrow = 0; }
            result.data[i] = (uint32_t)diff;
        }
        result.normalize();
        return result;
    }

    BigInt shiftLeft(int n) const {
        if (n == 0 || isZero()) return *this;
        BigInt result;
        result.data.clear();
        int wordShift = n / 32;
        int bitShift = n % 32;
        int newSize = data.size() + wordShift + (bitShift > 0 ? 1 : 0);
        result.data.resize(newSize, 0);

        if (bitShift == 0) {
            for (int i = 0; i < (int)data.size(); i++)
                result.data[i + wordShift] = data[i];
        } else {
            uint64_t carry = 0;
            for (int i = 0; i < (int)data.size(); i++) {
                uint64_t temp = ((uint64_t)data[i] << bitShift) | carry;
                result.data[i + wordShift] = (uint32_t)(temp & 0xFFFFFFFF);
                carry = temp >> 32;
            }
            if (carry && wordShift + data.size() < result.data.size())
                result.data[wordShift + data.size()] = (uint32_t)carry;
        }
        result.normalize();
        return result;
    }

    BigInt shiftRight(int n) const {
        if (n == 0 || isZero()) return *this;
        BigInt result;
        result.data.clear();
        int wordShift = n / 32;
        int bitShift = n % 32;
        if (wordShift >= (int)data.size()) return BigInt(0);

        int newSize = data.size() - wordShift;
        result.data.resize(newSize);
        for (int i = wordShift; i < (int)data.size(); i++)
            result.data[i - wordShift] = data[i];

        if (bitShift > 0) {
            for (int i = 0; i < (int)result.data.size(); i++) {
                result.data[i] = result.data[i] >> bitShift;
                if (i + 1 < (int)result.data.size())
                    result.data[i] |= (result.data[i + 1] & ((1U << bitShift) - 1)) << (32 - bitShift);
            }
        }
        result.normalize();
        return result;
    }

    BigInt operator*(const BigInt& other) const {
        if (isZero() || other.isZero()) return BigInt(0);
        BigInt result;
        result.data.assign(data.size() + other.data.size(), 0);

        for (int i = 0; i < (int)data.size(); i++) {
            uint64_t carry = 0;
            for (int j = 0; j < (int)other.data.size(); j++) {
                uint64_t prod = (uint64_t)data[i] * other.data[j];
                uint64_t sum = (uint64_t)result.data[i + j] + prod + carry;
                result.data[i + j] = (uint32_t)(sum & 0xFFFFFFFF);
                carry = sum >> 32;
            }
            if (carry) result.data[i + other.data.size()] += (uint32_t)carry;
        }
        result.normalize();
        return result;
    }

    void divMod(const BigInt& divisor, BigInt& quotient, BigInt& remainder) const {
        quotient = BigInt(0);
        remainder = BigInt(0);
        if (divisor.isZero()) return;
        if (*this < divisor) { remainder = *this; return; }
        if (divisor.isOne()) { quotient = *this; return; }

        if (divisor.data.size() == 1) {
            uint64_t div = divisor.data[0];
            uint64_t rem = 0;
            quotient.data.resize(data.size());
            for (int i = data.size() - 1; i >= 0; i--) {
                rem = (rem << 32) | data[i];
                quotient.data[i] = (uint32_t)(rem / div);
                rem = rem % div;
            }
            quotient.normalize();
            remainder.data.clear();
            remainder.data.push_back((uint32_t)rem);
            remainder.normalize();
            return;
        }

        int d = 0;
        uint32_t top = divisor.data.back();
        while ((top & 0x80000000) == 0) { top <<= 1; d++; }

        BigInt u = this->shiftLeft(d);
        BigInt v = divisor.shiftLeft(d);
        int n = v.data.size();
        int m = u.data.size() - n;
        quotient.data.resize(m + 1, 0);

        uint32_t vn_1 = v.data[n - 1];
        uint32_t vn_2 = v.data[n - 2];
        if ((int)u.data.size() <= m + n) u.data.resize(m + n + 1, 0);

        for (int j = m; j >= 0; j--) {
            uint64_t qhat, rhat;
            uint64_t u_jn = u.data[j + n];
            uint64_t u_jn_1 = u.data[j + n - 1];
            uint64_t u_jn_2 = u.data[j + n - 2];
            uint64_t dividend = (u_jn << 32) | u_jn_1;

            if (u_jn == vn_1) { qhat = 0xFFFFFFFF; rhat = u_jn_1 + vn_1; }
            else { qhat = dividend / vn_1; rhat = dividend % vn_1; }

            while (rhat < 0x100000000ULL && (qhat * vn_2 > (rhat << 32) + u_jn_2)) { qhat--; rhat += vn_1; }

            int64_t borrow = 0;
            for (int i = 0; i < n; i++) {
                uint64_t p = qhat * v.data[i];
                int64_t sub = (int64_t)u.data[j + i] - borrow - (p & 0xFFFFFFFF);
                u.data[j + i] = (uint32_t)sub;
                borrow = (p >> 32) - (sub >> 32);
            }

            int64_t sub = (int64_t)u.data[j + n] - borrow;
            u.data[j + n] = (uint32_t)sub;

            if (sub < 0) {
                qhat--;
                uint64_t carry = 0;
                for (int i = 0; i < n; i++) {
                    uint64_t sum = (uint64_t)u.data[j + i] + v.data[i] + carry;
                    u.data[j + i] = (uint32_t)sum;
                    carry = sum >> 32;
                }
                u.data[j + n] += (uint32_t)carry;
            }
            quotient.data[j] = (uint32_t)qhat;
        }
        quotient.normalize();
        remainder = u.shiftRight(d);
        remainder.normalize();
    }

    BigInt operator/(const BigInt& other) const { BigInt q, r; divMod(other, q, r); return q; }
    BigInt operator%(const BigInt& other) const { BigInt q, r; divMod(other, q, r); return r; }

    static BigInt mulMod(const BigInt& a, const BigInt& b, const BigInt& n) {
        if (n.isOne()) return BigInt(0);
        return (a * b) % n;
    }

    static BigInt powerMod(const BigInt& base, const BigInt& exp, const BigInt& n) {
        if (n.isOne()) return BigInt(0);
        BigInt result(1);
        BigInt b = base % n;
        int bits = exp.bitLength();
        for (int i = 0; i < bits; i++) {
            if (exp.getBit(i)) result = mulMod(result, b, n);
            b = mulMod(b, b, n);
        }
        return result;
    }

    static BigInt extendedGcd(const BigInt& a, const BigInt& b, BigInt& x, BigInt& y, bool& xNeg, bool& yNeg) {
        if (b.isZero()) { x = BigInt(1); y = BigInt(0); xNeg = false; yNeg = false; return a; }
        BigInt x1, y1; bool x1Neg, y1Neg;
        BigInt g = extendedGcd(b, a % b, x1, y1, x1Neg, y1Neg);
        x = y1; xNeg = y1Neg;
        BigInt q = a / b;
        BigInt qy1 = q * y1;
        if (x1Neg == y1Neg) {
            if (x1 >= qy1) { y = x1 - qy1; yNeg = x1Neg; }
            else { y = qy1 - x1; yNeg = !x1Neg; }
        } else { y = x1 + qy1; yNeg = x1Neg; }
        return g;
    }

    static BigInt modInverse(const BigInt& a, const BigInt& n) {
        BigInt x, y; bool xNeg, yNeg;
        BigInt g = extendedGcd(a % n, n, x, y, xNeg, yNeg);
        if (!g.isOne()) return BigInt(0);
        if (xNeg) x = n - (x % n);
        return x % n;
    }

    static BigInt gcd(const BigInt& a, const BigInt& b) {
        if (b.isZero()) return a;
        return gcd(b, a % b);
    }

    static BigInt randomBits(int bits, mt19937_64& rng) {
        BigInt result;
        result.data.clear();
        int numWords = (bits + 31) / 32;
        result.data.resize(numWords, 0);
        for (int i = 0; i < numWords; i++)
            result.data[i] = (uint32_t)(rng() & 0xFFFFFFFF);
        int extraBits = numWords * 32 - bits;
        if (extraBits > 0 && numWords > 0)
            result.data[numWords - 1] &= (0xFFFFFFFF >> extraBits);
        if (bits > 0) result.setBit(bits - 1);
        result.normalize();
        return result;
    }

    static BigInt randomRange(const BigInt& n, mt19937_64& rng) {
        if (n <= BigInt(2)) return BigInt(1);
        int bits = n.bitLength();
        BigInt result;
        do { result = randomBits(bits, rng); } while (result.isZero() || result >= n);
        return result;
    }

    friend ostream& operator<<(ostream& os, const BigInt& n);
};

ostream& operator<<(ostream& os, const BigInt& n) {
    if (n.isZero()) { os << "0"; return os; }
    int totalBits = n.bitLength();
    int totalHexDigits = (totalBits + 3) / 4;
    for (int i = 0; i < totalHexDigits; i++) {
        int bitPos = i * 4;
        int wordPos = bitPos / 32;
        int bitInWord = bitPos % 32;
        int digit = 0;
        if (wordPos < (int)n.data.size()) {
            digit = (n.data[wordPos] >> bitInWord) & 0xF;
            if (bitInWord > 28 && wordPos + 1 < (int)n.data.size())
                digit |= (n.data[wordPos + 1] & ((1 << (bitInWord - 28)) - 1)) << (32 - bitInWord);
        }
        os << (char)(digit < 10 ? '0' + digit : 'A' + digit - 10);
    }
    return os;
}

bool millerRabin(const BigInt& n, int k, mt19937_64& rng) {
    if (n < BigInt(2)) return false;
    if (n == BigInt(2)) return true;
    if (n.isEven()) return false;

    BigInt n_minus_1 = n - BigInt(1);
    BigInt d = n_minus_1;
    int r = 0;
    while (d.isEven()) { d = d.shiftRight(1); r++; }

    for (int i = 0; i < k; i++) {
        BigInt a = BigInt::randomRange(n_minus_1, rng);
        if (a < BigInt(2)) a = BigInt(2);
        BigInt x = BigInt::powerMod(a, d, n);
        if (x.isOne() || x == n_minus_1) continue;
        bool composite = true;
        for (int j = 0; j < r - 1; j++) {
            x = BigInt::mulMod(x, x, n);
            if (x == n_minus_1) { composite = false; break; }
        }
        if (composite) return false;
    }
    return true;
}

BigInt generatePrime(int bits, mt19937_64& rng) {
    while (true) {
        BigInt candidate = BigInt::randomBits(bits, rng);
        candidate.setBit(0);
        candidate.setBit(bits - 1);
        if (millerRabin(candidate, 25, rng)) return candidate;
    }
}

BigInt findPrimitiveRoot(const BigInt& p, mt19937_64& rng) {
    BigInt p_minus_1 = p - BigInt(1);
    vector<BigInt> factors;
    BigInt temp = p_minus_1;
    vector<uint64_t> smallPrimes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
    for (uint64_t sp : smallPrimes) {
        BigInt spBig(sp);
        if (temp % spBig == BigInt(0)) {
            factors.push_back(spBig);
            while (temp % spBig == BigInt(0)) temp = temp / spBig;
        }
    }
    if (temp > BigInt(1)) factors.push_back(temp);

    while (true) {
        BigInt g = BigInt::randomRange(p, rng);
        if (g <= BigInt(1)) continue;
        bool isPrimRoot = true;
        for (const BigInt& q : factors) {
            BigInt exp = p_minus_1 / q;
            if (BigInt::powerMod(g, exp, p).isOne()) { isPrimRoot = false; break; }
        }
        if (isPrimRoot) return g;
    }
}

bool verifyElGamalSignature(const BigInt& p, const BigInt& g, const BigInt& y,
                             const BigInt& m, const BigInt& r, const BigInt& h) {
    BigInt zero(0), one(1);
    BigInt p_minus_1 = p - one;
    if (r <= zero || r >= p) return false;
    if (h <= zero || h >= p_minus_1) return false;
    BigInt left = BigInt::powerMod(g, m, p);
    BigInt y_r = BigInt::powerMod(y, r, p);
    BigInt r_h = BigInt::powerMod(r, h, p);
    BigInt right = BigInt::mulMod(y_r, r_h, p);
    return left == right;
}

void generateValidSignature(const BigInt& p, const BigInt& g, const BigInt& x,
                            const BigInt& y, const BigInt& m,
                            BigInt& r, BigInt& h, mt19937_64& rng) {
    BigInt p_minus_1 = p - BigInt(1);
    BigInt k;
    do { k = BigInt::randomRange(p_minus_1, rng); } while (BigInt::gcd(k, p_minus_1) != BigInt(1));
    r = BigInt::powerMod(g, k, p);
    BigInt k_inv = BigInt::modInverse(k, p_minus_1);
    BigInt xr = BigInt::mulMod(x, r, p_minus_1);
    BigInt m_mod = m % p_minus_1;
    BigInt diff = (m_mod >= xr) ? m_mod - xr : p_minus_1 - (xr - m_mod);
    h = BigInt::mulMod(diff, k_inv, p_minus_1);
    if (h.isZero()) h = BigInt(1);
}

void writeTestCase(const string& filename, const BigInt& p, const BigInt& g,
                   const BigInt& y, const BigInt& m, const BigInt& r, const BigInt& h) {
    ofstream out(filename);
    out << p << "\n" << g << "\n" << y << "\n" << m << "\n" << r << "\n" << h << "\n";
    out.close();
}

string readFile(const string& filename) {
    ifstream file(filename);
    if (!file) return "";
    string content;
    getline(file, content);
    return content;
}

int main(int argc, char* argv[]) {
    string outputDir = ".";
    int numTests = 200;  // Generate 200 tests by default

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if ((arg == "--dir" || arg == "-d") && i + 1 < argc) outputDir = argv[++i];
        else if ((arg == "--num" || arg == "-n") && i + 1 < argc) numTests = atoi(argv[++i]);
        else if (arg == "--help" || arg == "-h") {
            cout << "Usage: " << argv[0] << " [options]\n";
            cout << "Options:\n";
            cout << "  -d, --dir DIR    Output directory (default: .)\n";
            cout << "  -n, --num N      Number of tests (default: 200)\n";
            cout << "  -h, --help       Show this help\n";
            return 0;
        }
    }

    mt19937_64 rng(time(nullptr));

    cout << "=== GENERATING " << numTests << " TEST CASES ===" << endl;

    // More complex distribution with larger bit ranges
    // Section 4.4: 256-512: 10%, 128-256: 20%, 64-128: 30%, <=64: 40%
    struct BitRange {
        int minBits, maxBits, percent;
        string name;
    };

    vector<BitRange> ranges = {
        {384, 512, 5, "Very Large (384-512 bits)"},
        {256, 384, 5, "Large (256-384 bits)"},
        {192, 256, 10, "Medium-Large (192-256 bits)"},
        {128, 192, 10, "Medium (128-192 bits)"},
        {96, 128, 15, "Medium-Small (96-128 bits)"},
        {64, 96, 15, "Small-Medium (64-96 bits)"},
        {32, 64, 20, "Small (32-64 bits)"},
        {16, 32, 10, "Very Small (16-32 bits)"},
        {8, 16, 10, "Tiny (8-16 bits)"}
    };

    int testNum = 0;
    int totalValid = 0, totalInvalid = 0;

    for (const auto& range : ranges) {
        int count = (numTests * range.percent) / 100;
        if (count == 0) count = 1;

        cout << "\n" << range.name << " (" << count << " tests):" << endl;

        for (int i = 0; i < count && testNum < numTests; i++) {
            int bits = range.minBits + (rng() % (range.maxBits - range.minBits + 1));
            if (bits < 8) bits = 8;

            cout << "  Test " << setw(3) << testNum << " (" << setw(3) << bits << " bits)..." << flush;

            BigInt p = generatePrime(bits, rng);
            BigInt p_minus_1 = p - BigInt(1);
            BigInt g = findPrimitiveRoot(p, rng);
            BigInt x = BigInt::randomRange(p_minus_1, rng);
            if (x.isZero()) x = BigInt(1);
            BigInt y = BigInt::powerMod(g, x, p);
            BigInt m = BigInt::randomRange(p_minus_1, rng);
            if (m.isZero()) m = BigInt(1);

            BigInt r, h;

            // Alternate valid/invalid like professor's tests
            bool valid = (testNum % 2 == 0);

            int attempts = 0;
            while (attempts < 10) {
                if (valid) {
                    generateValidSignature(p, g, x, y, m, r, h, rng);
                    if (verifyElGamalSignature(p, g, y, m, r, h)) break;
                } else {
                    // Generate invalid signature
                    r = BigInt::randomRange(p, rng);
                    if (r.isZero()) r = BigInt(1);
                    h = BigInt::randomRange(p_minus_1, rng);
                    if (h.isZero()) h = BigInt(1);
                    if (!verifyElGamalSignature(p, g, y, m, r, h)) break;
                }
                attempts++;
            }

            if (attempts >= 10) {
                cout << " (retry)..." << flush;
                i--;
                continue;
            }

            string inputFile = outputDir + "/test_" + to_string(testNum) + ".inp";
            writeTestCase(inputFile, p, g, y, m, r, h);

            string outputFile = outputDir + "/test_" + to_string(testNum) + ".out";
            ofstream out(outputFile);
            out << (valid ? "1" : "0") << endl;
            out.close();

            if (valid) totalValid++; else totalInvalid++;
            cout << " " << (valid ? "VALID" : "INVALID") << endl;
            testNum++;
        }
    }

    cout << "\n=== GENERATION COMPLETE ===" << endl;
    cout << "Total tests: " << testNum << endl;
    cout << "Valid signatures: " << totalValid << endl;
    cout << "Invalid signatures: " << totalInvalid << endl;
    cout << "Output directory: " << outputDir << "/" << endl;

    // Auto-verify
    cout << "\n=== VERIFYING TEST CASES ===" << endl;
    cout << "Compiling main.cpp..." << endl;
    int compileResult = system("g++ -std=c++17 -O3 -o main_verify main.cpp 2>/dev/null");

    if (compileResult != 0) {
        cerr << "Failed to compile main.cpp, skipping verification." << endl;
        return 0;
    }

    int passed = 0, failed = 0;
    for (int i = 0; i < testNum; i++) {
        string inputFile = outputDir + "/test_" + to_string(i) + ".inp";
        string expectedFile = outputDir + "/test_" + to_string(i) + ".out";
        string actualFile = outputDir + "/test_actual_" + to_string(i) + ".out";

        string cmd = "./main_verify " + inputFile + " " + actualFile + " 2>/dev/null";
        system(cmd.c_str());

        string expected = readFile(expectedFile);
        string actual = readFile(actualFile);

        // Trim whitespace
        while (!expected.empty() && (expected.back() == '\n' || expected.back() == '\r')) expected.pop_back();
        while (!actual.empty() && (actual.back() == '\n' || actual.back() == '\r')) actual.pop_back();

        if (expected == actual) {
            passed++;
        } else {
            cout << "Test " << i << ": FAIL (expected=" << expected << ", got=" << actual << ")" << endl;
            failed++;
        }
        remove(actualFile.c_str());
    }

    cout << "\n=== VERIFICATION SUMMARY ===" << endl;
    cout << "Passed: " << passed << "/" << testNum << endl;
    cout << "Failed: " << failed << "/" << testNum << endl;
    if (failed == 0) cout << "All tests passed!" << endl;

    return 0;
}
