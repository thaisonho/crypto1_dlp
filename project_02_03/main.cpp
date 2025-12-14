#pragma GCC optimize("O3")
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// Big integer for non-negative numbers.
// input format h_0*16^0 + h_1*16^1 + ... (first char is LSB)
class BigInt {
private:
    vector<uint32_t> data;

    void normalize() {
        while (data.size() > 1 && data.back() == 0) data.pop_back();
        if (data.empty()) data.push_back(0);
    }

public:
    BigInt() { data.push_back(0); }

    BigInt(uint64_t val) {
        data.push_back((uint32_t)(val & 0xFFFFFFFFu));
        if (val > 0xFFFFFFFFu) data.push_back((uint32_t)(val >> 32));
        normalize();
    }

    // input format h_0*16^0 + h_1*16^1 + ... (first char is LSB)
    BigInt(const string& hex) {
        if (hex.empty() || hex == "0") {
            data.push_back(0);
            return;
        }

        // Estimate size to avoid reallocations
        // Each hex digit is 4 bits. 32 bits per word.
        // So roughly hex.length() / 8 words.
        data.reserve((hex.length() + 7) / 8);
        data.push_back(0); // Start with 0

        for (int i = 0; i < (int)hex.length(); i++) {
            char c = hex[i];  // h_i is at position i (left to right)
            int digit;
            if (c >= '0' && c <= '9') digit = c - '0';
            else if (c >= 'A' && c <= 'F') digit = c - 'A' + 10;
            else if (c >= 'a' && c <= 'f') digit = c - 'a' + 10;
            else continue;

            // h_i contributes to bit positions [4*i, 4*i+3]
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

    bool getBit(int pos) const {
        int wordPos = pos / 32;
        int bitPos = pos % 32;
        if (wordPos >= (int)data.size()) return false;
        return (data[wordPos] >> bitPos) & 1u;
    }

    int bitLength() const {
        if (isZero()) return 0;
        int len = (int)(data.size() - 1) * 32;
        uint32_t top = data.back();
        while (top) {
            len++;
            top >>= 1;
        }
        return len;
    }

    bool operator==(const BigInt& other) const { return data == other.data; }
    bool operator!=(const BigInt& other) const { return !(*this == other); }

    bool operator<(const BigInt& other) const {
        if (data.size() != other.data.size()) return data.size() < other.data.size();
        for (int i = (int)data.size() - 1; i >= 0; i--) {
            if (data[i] != other.data[i]) return data[i] < other.data[i];
        }
        return false;
    }
    bool operator>=(const BigInt& other) const { return !(*this < other); }

    BigInt operator+(const BigInt& other) const {
        BigInt result;
        result.data.clear(); // Clear initial 0
        uint64_t carry = 0;
        int maxSize = max((int)data.size(), (int)other.data.size());
        result.data.reserve(maxSize + 1);
        for (int i = 0; i < maxSize || carry; i++) {
            uint64_t sum = carry;
            if (i < (int)data.size()) sum += data[i];
            if (i < (int)other.data.size()) sum += other.data[i];
            result.data.push_back((uint32_t)(sum & 0xFFFFFFFFu));
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
            if (diff < 0) {
                diff += 0x100000000LL;
                borrow = 1;
            } else {
                borrow = 0;
            }
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
        int newSize = (int)data.size() + wordShift + (bitShift ? 1 : 0);
        result.data.resize(newSize, 0);
        if (bitShift == 0) {
            for (int i = 0; i < (int)data.size(); i++) result.data[i + wordShift] = data[i];
        } else {
            uint64_t carry = 0;
            for (int i = 0; i < (int)data.size(); i++) {
                uint64_t temp = ((uint64_t)data[i] << bitShift) | carry;
                result.data[i + wordShift] = (uint32_t)(temp & 0xFFFFFFFFu);
                carry = temp >> 32;
            }
            if (carry) {
                if (wordShift + (int)data.size() < (int)result.data.size()) {
                    result.data[wordShift + (int)data.size()] = (uint32_t)carry;
                } else {
                    result.data.push_back((uint32_t)carry);
                }
            }
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
        int newSize = (int)data.size() - wordShift;
        result.data.resize(newSize);
        for (int i = wordShift; i < (int)data.size(); i++) result.data[i - wordShift] = data[i];
        if (bitShift > 0) {
            for (int i = 0; i < (int)result.data.size(); i++) {
                result.data[i] >>= bitShift;
                if (i + 1 < (int)result.data.size()) {
                    result.data[i] |= (result.data[i + 1] & ((1u << bitShift) - 1u)) << (32 - bitShift);
                }
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
                uint64_t prod = (uint64_t)data[i] * (uint64_t)other.data[j];
                uint64_t sum = (uint64_t)result.data[i + j] + prod + carry;
                result.data[i + j] = (uint32_t)(sum & 0xFFFFFFFFu);
                carry = sum >> 32;
            }
            if (carry) result.data[i + (int)other.data.size()] += (uint32_t)carry;
        }
        result.normalize();
        return result;
    }

    void divMod(const BigInt& divisor, BigInt& quotient, BigInt& remainder) const {
        quotient = BigInt(0);
        remainder = BigInt(0);
        if (divisor.isZero()) return;
        if (*this < divisor) {
            remainder = *this;
            return;
        }
        if (divisor.isOne()) {
            quotient = *this;
            return;
        }

        // Fast path: single-word divisor
        if (divisor.data.size() == 1) {
            uint64_t div = divisor.data[0];
            uint64_t rem = 0;
            quotient.data.resize(data.size());
            for (int i = (int)data.size() - 1; i >= 0; i--) {
                rem = (rem << 32) | data[i];
                quotient.data[i] = (uint32_t)(rem / div);
                rem %= div;
            }
            quotient.normalize();
            remainder.data.clear();
            remainder.data.push_back((uint32_t)rem);
            remainder.normalize();
            return;
        }

        // Knuth's Algorithm D
        // Normalize
        int d = 0;
        uint32_t top = divisor.data.back();
        while ((top & 0x80000000u) == 0) {
            top <<= 1;
            d++;
        }

        BigInt u = this->shiftLeft(d);
        BigInt v = divisor.shiftLeft(d);

        int n = (int)v.data.size();
        int m = (int)u.data.size() - n;
        quotient.data.resize(m + 1, 0);

        uint32_t vn_1 = v.data[n - 1];
        uint32_t vn_2 = v.data[n - 2];

        if ((int)u.data.size() <= m + n) u.data.resize(m + n + 1, 0);

        for (int j = m; j >= 0; j--) {
            uint64_t qhat;
            uint64_t rhat;
            uint64_t u_jn = u.data[j + n];
            uint64_t u_jn_1 = u.data[j + n - 1];
            uint64_t u_jn_2 = u.data[j + n - 2];

            uint64_t dividend = (u_jn << 32) | u_jn_1;

            if (u_jn == vn_1) {
                qhat = 0xFFFFFFFFu;
                rhat = u_jn_1 + vn_1;
            } else {
                qhat = dividend / vn_1;
                rhat = dividend % vn_1;
            }

            while (rhat < 0x100000000ULL && (qhat * vn_2 > (rhat << 32) + u_jn_2)) {
                qhat--;
                rhat += vn_1;
            }

            // Multiply and subtract
            int64_t borrow = 0;
            for (int i = 0; i < n; i++) {
                uint64_t p = qhat * (uint64_t)v.data[i];
                int64_t sub = (int64_t)u.data[j + i] - borrow - (int64_t)(p & 0xFFFFFFFFu);
                u.data[j + i] = (uint32_t)sub;
                borrow = (int64_t)(p >> 32) - (sub >> 32);
            }

            // Handle the top digit of u (u_{j+n})
            int64_t subTop = (int64_t)u.data[j + n] - borrow;
            u.data[j + n] = (uint32_t)subTop;

            // If we subtracted too much, add back
            if (subTop < 0) {
                qhat--;
                uint64_t carry = 0;
                for (int i = 0; i < n; i++) {
                    uint64_t sum = (uint64_t)u.data[j + i] + (uint64_t)v.data[i] + carry;
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

    // Modular multiplication using optimized multiplication + reduction

    BigInt operator/(const BigInt& other) const {
        BigInt q, r;
        divMod(other, q, r);
        return q;
    }

    BigInt operator%(const BigInt& other) const {
        BigInt q, r;
        divMod(other, q, r);
        return r;
    }

    static BigInt mulMod(const BigInt& a, const BigInt& b, const BigInt& n) {
        if (n.isOne()) return BigInt(0);

        // Use optimized multiplication
        BigInt product = a * b;

        // Fast reduction for result
        return product % n;
    }

    static BigInt powerMod(const BigInt& base, const BigInt& exp, const BigInt& n) {
        if (n.isOne()) return BigInt(0);

        BigInt result(1);
        BigInt b = base % n;

        int bits = exp.bitLength();
        for (int i = 0; i < bits; i++) {
            if (exp.getBit(i)) {
                result = mulMod(result, b, n);
            }
            b = mulMod(b, b, n);
        }

        return result;
    }

    friend ostream& operator<<(ostream& os, const BigInt& n);
};

ostream& operator<<(ostream& os, const BigInt& n) {
    if (n.isZero()) {
        os << "0";
        return os;
    }

    // Calculate total number of hex digits needed
    int totalBits = n.bitLength();
    int totalHexDigits = (totalBits + 3) / 4;  // Round up to nearest hex digit

    // Output each hex digit from LSB to MSB (left to right)
    for (int i = 0; i < totalHexDigits; i++) {
        // Extract 4 bits starting at position i*4
        int bitPos = i * 4;
        int wordPos = bitPos / 32;
        int bitInWord = bitPos % 32;

        int digit = 0;
        if (wordPos < (int)n.data.size()) {
            digit = (n.data[wordPos] >> bitInWord) & 0xF;

            // Handle case where hex digit spans two words
            if (bitInWord > 28 && wordPos + 1 < (int)n.data.size()) {
                int bitsFromNextWord = bitInWord - 28;
                digit |= (n.data[wordPos + 1] & ((1 << bitsFromNextWord) - 1)) << (32 - bitInWord);
            }
        }

        // Output as hex character
        if (digit < 10) {
            os << (char)('0' + digit);
        }
        else {
            os << (char)('A' + digit - 10);
        }
    }

    return os;
}

// ElGamal decryption
// Input (5 lines): p, g, x, c1, c2   (all in hex little-endian-by-nibble)
// Output (1 line): m
//
// Decryption formula:
//   s   = c1^x mod p
//   inv = s^(p-2) mod p  (Fermat's little theorem, since p is prime)
//   m   = c2 * inv mod p
int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>\n";
        return 1;
    }

    ifstream in(argv[1]);
    if (!in) {
        cerr << "Cannot open input file: " << argv[1] << "\n";
        return 1;
    }

    string pHex, gHex, xHex, c1Hex, c2Hex;
    in >> pHex >> gHex >> xHex >> c1Hex >> c2Hex;
    in.close();

    BigInt p(pHex);
    BigInt g(gHex);
    (void)g; // g is not needed for decryption computation (kept for input format)
    BigInt x(xHex);
    BigInt c1(c1Hex);
    BigInt c2(c2Hex);

    // s = c1^x mod p
    BigInt s = BigInt::powerMod(c1, x, p);
    // inv = s^(p-2) mod p
    BigInt pMinus2 = p - BigInt(2);
    BigInt inv = BigInt::powerMod(s, pMinus2, p);
    // m = c2 * inv mod p
    BigInt m = BigInt::mulMod(c2 % p, inv, p);

    ofstream out(argv[2]);
    if (!out) {
        cerr << "Cannot open output file: " << argv[2] << "\n";
        return 1;
    }
    out << m << "\n";
    out.close();
    return 0;
}


