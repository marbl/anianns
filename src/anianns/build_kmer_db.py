import gzip
import struct

def load_kmer_sets_shared_k(path: str) -> tuple[int, dict[str, set[int]]]:
    try:
        with gzip.open(path, "rb") as f:
            # Read shared k
            raw_k = f.read(1)
            if len(raw_k) != 1:
                raise ValueError("Failed to read k-mer length (1 byte expected)")

            k = struct.unpack("<B", raw_k)[0]

            # Read number of entries
            raw_num = f.read(2)
            if len(raw_num) != 2:
                raise ValueError("Failed to read number of entries (2 bytes expected)")

            num_sets = struct.unpack("<H", raw_num)[0]
            result = {}

            for i in range(num_sets):
                # Name length
                raw_len = f.read(2)
                if len(raw_len) != 2:
                    raise ValueError(f"[Entry {i}] Failed to read name length")

                name_len = struct.unpack("<H", raw_len)[0]

                # Name string
                name_bytes = f.read(name_len)
                if len(name_bytes) != name_len:
                    raise ValueError(f"[Entry {i}] Incomplete name data")

                name = name_bytes.decode("utf-8")

                # Number of kmers
                raw_kmer_count = f.read(4)
                if len(raw_kmer_count) != 4:
                    raise ValueError(f"[{name}] Failed to read number of kmers")

                num_kmers = struct.unpack("<I", raw_kmer_count)[0]

                # Kmers
                kmers = set()
                for j in range(num_kmers):
                    raw_kmer = f.read(4)
                    if len(raw_kmer) != 4:
                        raise ValueError(
                            f"[{name}] Failed to read kmer {j+1}/{num_kmers}"
                        )
                    kmers.add(struct.unpack("<i", raw_kmer)[0])

                result[name] = kmers

        return k, result

    except (OSError, struct.error, ValueError, UnicodeDecodeError) as e:
        print(f"[ERROR] Failed to load file '{path}': {e}")
        return None, {}


def save_kmer_sets_shared_k(
    kmer_dict: dict[str, set[int]], k: int, output_path: str
) -> None:
    with gzip.open(output_path, "wb") as f:
        f.write(struct.pack("<B", k))  # uint8: shared k-mer length
        f.write(struct.pack("<H", len(kmer_dict)))  # uint16: number of entries

        for name, kmers in kmer_dict.items():
            name_bytes = name.encode("utf-8")
            f.write(struct.pack("<H", len(name_bytes)))  # uint16: name length
            f.write(name_bytes)  # name string
            f.write(struct.pack("<I", len(kmers)))  # uint32: number of kmers

            for kmer in kmers:
                f.write(struct.pack("<i", kmer))  # int32: each kmer
