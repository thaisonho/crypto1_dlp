#!/usr/bin/env python3
"""
Test generator and validator for Project 02_03: ElGamal (decrypt).

Input format (5 lines, hex little-endian-by-nibble):
  p
  g
  x
  c1
  c2

Output format (1 line, hex little-endian-by-nibble):
  m

We generate valid ElGamal ciphertexts from random (p, g, x, k, m):
  y  = g^x mod p
  c1 = g^k mod p
  c2 = m * y^k mod p
and expected plaintext is m.
"""

import os
import subprocess
import tempfile
import random

import gmpy2
from gmpy2 import mpz, next_prime, powmod


def int_to_lsb_hex(n: int) -> str:
    """Integer -> hex string in project format (LSB hex digit first)."""
    if n == 0:
        return "0"
    out = []
    n = int(n)
    while n > 0:
        digit = n & 0xF
        if digit < 10:
            out.append(chr(ord("0") + digit))
        else:
            out.append(chr(ord("A") + digit - 10))
        n >>= 4
    return "".join(out)


def lsb_hex_to_int(hex_str: str) -> int:
    """Project format hex -> integer."""
    result = 0
    for i, c in enumerate(hex_str.strip()):
        if "0" <= c <= "9":
            digit = ord(c) - ord("0")
        elif "A" <= c <= "F":
            digit = ord(c) - ord("A") + 10
        elif "a" <= c <= "f":
            digit = ord(c) - ord("a") + 10
        else:
            continue
        result += digit * (16**i)
    return result


def random_prime(bits: int) -> mpz:
    """Generate a prime of ~bits bits using next_prime."""
    if bits < 8:
        bits = 8
    # ensure top bit set so it is truly "bits" sized most of the time
    n = (mpz(1) << (bits - 1)) | mpz(random.getrandbits(bits - 1)) | 1
    p = next_prime(n)
    # If next_prime jumped bit-length, try again a few times.
    for _ in range(10):
        if p.bit_length() == bits:
            return p
        n = (mpz(1) << (bits - 1)) | mpz(random.getrandbits(bits - 1)) | 1
        p = next_prime(n)
    return p


def generate_case(bits: int):
    """
    Returns (p, g, x, c1, c2, m_expected) as mpz/int.
    Constraints: g, x, c1, c2 < p and 0 <= m < p.
    """
    p = random_prime(bits)
    p_int = int(p)

    g = mpz(random.randrange(2, p_int - 1))
    x = mpz(random.randrange(1, p_int - 1))
    k = mpz(random.randrange(1, p_int - 1))
    m = mpz(random.randrange(0, p_int))

    y = powmod(g, x, p)
    c1 = powmod(g, k, p)
    c2 = (m * powmod(y, k, p)) % p

    return p, g, x, c1, c2, m


def create_input(p: mpz, g: mpz, x: mpz, c1: mpz, c2: mpz) -> str:
    return (
        int_to_lsb_hex(int(p)) + "\n"
        + int_to_lsb_hex(int(g)) + "\n"
        + int_to_lsb_hex(int(x)) + "\n"
        + int_to_lsb_hex(int(c1)) + "\n"
        + int_to_lsb_hex(int(c2)) + "\n"
    )


def run_cpp_program(input_content: str, executable: str) -> str | None:
    inp_fd, inp_path = tempfile.mkstemp(suffix=".inp", text=True)
    out_fd, out_path = tempfile.mkstemp(suffix=".out", text=True)
    os.close(inp_fd)
    os.close(out_fd)

    try:
        with open(inp_path, "w") as f:
            f.write(input_content)

        result = subprocess.run(
            [executable, inp_path, out_path],
            capture_output=True,
            text=True,
            timeout=60,
        )
        if result.returncode != 0:
            print(f"C++ program error: {result.stderr}")
            return None

        with open(out_path, "r") as f:
            return f.read().strip()
    except subprocess.TimeoutExpired:
        print("C++ program timed out")
        return None
    finally:
        try:
            os.unlink(inp_path)
            os.unlink(out_path)
        except OSError:
            pass


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Test generator for ElGamal decryption (Project 02_03)")
    parser.add_argument("--num-tests", type=int, default=100, help="Number of tests to generate")
    parser.add_argument("--executable", type=str, default="./main", help="Path to C++ executable")
    parser.add_argument("--save-tests", type=str, default=None, help="Directory to save generated tests (.inp/.out)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--max-bits", type=int, default=512, help="Maximum bit size for primes")
    args = parser.parse_args()

    if not os.path.exists(args.executable):
        print(f"Executable not found: {args.executable}")
        print("Please compile first: g++ -O3 -std=c++17 -o main main.cpp")
        return 1

    if args.save_tests:
        os.makedirs(args.save_tests, exist_ok=True)

    # Ensure we always have a valid bit-range for generation logic.
    # (Our prime generator uses at least 8 bits anyway.)
    max_bits = max(8, int(args.max_bits))

    # Match the rubric distribution:
    # 40%: <=64, 30%: (64,128], 20%: (128,256], 10%: (256,512]
    bit_sizes: list[int] = []

    # Desired counts per bucket (we'll reassign buckets that are not feasible).
    n1 = args.num_tests * 40 // 100
    n2 = args.num_tests * 30 // 100
    n3 = args.num_tests * 20 // 100
    n4 = args.num_tests * 10 // 100

    # Bucket definitions: (low, high, count)
    # We will generate from a bucket only if low <= high after clamping to max_bits.
    buckets = [
        (8,  min(64,  max_bits), n1),
        (65, min(128, max_bits), n2),
        (129, min(256, max_bits), n3),
        (257, min(512, max_bits), n4),
    ]

    # Approach B: validate each bucket range; if invalid, fall back to the last valid range.
    last_valid = (8, min(64, max_bits))
    for lo, hi, cnt in buckets:
        if cnt <= 0:
            continue
        if lo <= hi:
            last_valid = (lo, hi)
            bit_sizes.extend([random.randint(lo, hi) for _ in range(cnt)])
        else:
            # Fallback to the nearest earlier feasible bucket range.
            fb_lo, fb_hi = last_valid
            bit_sizes.extend([random.randint(fb_lo, fb_hi) for _ in range(cnt)])

    while len(bit_sizes) < args.num_tests:
        bit_sizes.append(random.randint(8, min(64, max_bits)))
    bit_sizes = bit_sizes[: args.num_tests]
    random.shuffle(bit_sizes)

    passed = 0
    failed = 0
    errors = 0

    print(f"Running {args.num_tests} tests...")
    print("-" * 60)

    for i, bits in enumerate(bit_sizes):
        try:
            p, g, x, c1, c2, expected_m = generate_case(bits)
            input_content = create_input(p, g, x, c1, c2)

            if args.save_tests:
                inp_path = os.path.join(args.save_tests, f"test_{i:03d}.inp")
                out_path = os.path.join(args.save_tests, f"test_{i:03d}.out")
                with open(inp_path, "w") as f:
                    f.write(input_content)
                with open(out_path, "w") as f:
                    f.write(int_to_lsb_hex(int(expected_m)) + "\n")

            cpp_out = run_cpp_program(input_content, args.executable)
            if cpp_out is None:
                errors += 1
                print(f"Test {i:3d}: ERROR (bits={bits})")
                continue

            expected_hex = int_to_lsb_hex(int(expected_m))
            if cpp_out.strip() == expected_hex:
                passed += 1
                if args.verbose:
                    print(f"Test {i:3d}: PASS (bits={bits})")
            else:
                failed += 1
                print(f"Test {i:3d}: FAIL (bits={bits})")
                print(f"  Expected: {expected_hex}")
                print(f"  Got     : {cpp_out.strip()}")
                if args.verbose:
                    print(f"  Input:\n{input_content}")
        except Exception as e:
            errors += 1
            print(f"Test {i:3d}: EXCEPTION - {e}")
            import traceback

            traceback.print_exc()

    print("-" * 60)
    print(f"Results: {passed} passed, {failed} failed, {errors} errors out of {args.num_tests} tests")
    if failed == 0 and errors == 0:
        print("All tests passed!")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())


