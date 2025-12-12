## Usage

### Running the Main Program

```bash
./main <input_file> <output_file>
```

**Example:**
```bash
./main test.inp test.out
```

### Input Format (`test.inp`)

The input file contains 6 lines, each with a non-negative integer in **little-endian hexadecimal** format:

```
<p>     Line 1: Prime number p
<g>     Line 2: Generator g (primitive root modulo p)
<y>     Line 3: Public key y
<m>     Line 4: Message m (where m < p - 1)
<r>     Line 5: Signature component r
<h>     Line 6: Signature component h
```

**Note**: All values are written as lowercase hexadecimal digits in little-endian format.

**Example** (test_00.inp):
```
16
71
F5
24
A5
A5
```

### Output Format (`test.out`)

Single line containing:
- `1` if the signature is **valid**
- `0` if the signature is **invalid**

## Verification Algorithm

Given public key `v = (p, g, y)`, message `m`, and signature `c = (r, h)`:

1. **Range Check**: Verify `0 < r < p` and `0 < h < p - 1`
2. **Signature Check**: Verify `g^m ≡ y^r × r^h (mod p)`

If both conditions are satisfied, output `1` (valid). Otherwise, output `0` (invalid).

## Test Generator

### Generate Test Cases

```bash
# Generate 200 test cases (default)
./test_generator

# Generate custom number of tests
./test_generator --num 100

# Generate tests in specific directory
./test_generator --dir ./my_tests

# Show help
./test_generator --help
```

## Generated Files

For each test case `N`, the generator creates:
- `test_N.inp` - Input file
- `test_N.out` - Expected output file

## Comparing with Sample Test Cases

### Run Against Professor's Sample Tests

```bash
# Test against the 10 sample tests with expected outputs
echo "=== Testing against sample tests ==="
for i in $(seq 0 9); do
    testnum=$(printf "%02d" $i)
    ./main ../testcases/project_02_04/test_${testnum}.inp my_output.out
    expected=$(cat ../testcases/project_02_04/test_${testnum}.out | tr -d '\n\r')
    actual=$(cat my_output.out | tr -d '\n\r')
    if [ "$expected" = "$actual" ]; then
        echo "Test $testnum: PASS (expected=$expected, got=$actual)"
    else
        echo "Test $testnum: FAIL (expected=$expected, got=$actual)"
    fi
done
```

### Run All 20 Sample Input Tests

```bash
# Run tests 00-19 (10-19 don't have expected outputs)
for i in $(seq 0 19); do
    testnum=$(printf "%02d" $i)
    ./main ../testcases/project_02_04/test_${testnum}.inp output_${testnum}.out
    echo "Test $testnum: $(cat output_${testnum}.out)"
done
```

## Sample Test Cases

### Test 01 (Valid Signature)

**Input** (`test_00.inp`):
```
16      (p = 0x61 = 97)
71      (g = 0x17 = 23)
F5      (y = 0x5F = 95)
24      (m = 0x42 = 66)
A5      (r = 0x5A = 90)
A5      (h = 0x5A = 90)
```

**Expected Output**: `1` (valid)

**Verification**:
- g^m = 23^66 = 85 (mod 97)
- y^r × r^h = 95^90 × 90^90 = 47 × 8 = 85 (mod 97)
- Since 85 = 85, signature is valid

### Test 02 (Invalid Signature)

**Input** (`test_01.inp`):
```
35      (p = 0x53 = 83)
31      (g = 0x13 = 19)
84      (y = 0x48 = 72)
24      (m = 0x42 = 66)
93      (r = 0x39 = 57)
E4      (h = 0x4E = 78)
```

**Expected Output**: `0` (invalid)

**Verification**:
- g^m = 19^66 = 78 (mod 83)
- y^r × r^h = 66 × 49 = 80 (mod 83)
- Since 78 ≠ 80, signature is invalid

## Quick Start

```bash
cd project_02_04

# 1. Compile
g++ -std=c++17 -O3 -o main main.cpp
g++ -std=c++17 -O3 -o test_generator test_generator.cpp

# 2. Generate and verify custom tests
./test_generator --num 100
```

## Author

Implementation for Cryptography Course - Exercise 4: ElGamal Digital Signature System
