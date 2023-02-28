import collections
from functools import cache
from collections.abc import Iterator
from datetime import datetime
from pathlib import Path
from typing import Literal, TypeAlias
from time import perf_counter

from libs.slide_min import Index, SlideMin, Value

Minimizer: TypeAlias = int
RefId: TypeAlias = str
ReadId: TypeAlias = str
Strand = Literal[0, 1]

CURRENT_FOLDER_DIR = Path(__file__).resolve().parent

NUCLEOTIDE_TO_INTEGER_MAPPING: dict[
    Literal["A", "C", "G", "T"], Literal[0, 1, 2, 3]
] = {"A": 0, "C": 1, "G": 2, "T": 3}


class Minimap:
    read_seq_file_path: Path
    ref_seq_file_path: Path
    mapping_result_output_file_path: Path
    window_size: int
    k_mer_size: int
    k_mer_hash_mod: int

    def __init__(
        self,
        read_seq_file_path: Path,
        ref_seq_file_path: Path,
        window_size: int,
        k_mer_size: int,
        mapping_result_output_file_path: Path | None = None,
    ) -> None:
        self.read_seq_file_path = read_seq_file_path
        self.ref_seq_file_path = ref_seq_file_path
        self.mapping_result_output_file_path = (
            CURRENT_FOLDER_DIR
            / f'mapping_result_{datetime.now().isoformat(timespec="seconds")}.tsv'
            if mapping_result_output_file_path is None
            else mapping_result_output_file_path
        )
        self.k_mer_size = k_mer_size
        self.k_mer_hash_mod = 1 << ((k_mer_size - 1) * 2)  # = 4 ** (k_mer_size - 1)
        self.window_size = window_size

    def parse_ref_seq_file(self) -> dict[RefId, str]:
        ref_seq_dict: dict[RefId, str] = {}
        with self.ref_seq_file_path.open("r") as f:
            tmp_ref_id: RefId | None = None
            seq = ""

            while True:
                line = f.readline().rstrip()
                if len(line) == 0:
                    break
                if line.startswith(">"):
                    if seq != "":
                        ref_seq_dict[tmp_ref_id] = seq
                        seq = ""
                    tmp_ref_id = line.split()[0][1:]
                else:
                    seq += line
            ref_seq_dict[tmp_ref_id] = seq
        return ref_seq_dict

    def parse_read_seq_file(self) -> dict[ReadId, str]:
        read_seq_dict: dict[ReadId, str] = {}
        with self.read_seq_file_path.open("r") as f:
            while True:
                read_id = f.readline().rstrip()[1:]
                read_array = f.readline().rstrip()
                f.readline()
                f.readline()
                if len(read_id) == 0:
                    break
                read_seq_dict[read_id] = read_array
        return read_seq_dict

    def sketch_minimizers(self, seq: str) -> Iterator[tuple[Value, Index, Strand]]:
        window_size = self.window_size
        k_mer_size = self.k_mer_size
        minimizers_set: set[tuple[Index, Strand]] = set()
        seq_length = len(seq)
        k_mer_init = seq[:k_mer_size]
        k_mer_hashes_strand_0 = SlideMin(self.seq2hash(k_mer_init), window_size)
        k_mer_hashes_strand_1 = SlideMin(
            self.seq2hash(k_mer_init, strand=1), window_size
        )
        for i in range(window_size - 1):
            last_k_mer_hash_strand_0_value = k_mer_hashes_strand_0.last_value
            last_k_mer_hash_strand_1_value = k_mer_hashes_strand_1.last_value
            (
                new_k_mer_hash_strand_0_value,
                new_k_mer_hash_strand_1_value,
            ) = self.calc_hash(
                last_k_mer_hash_strand_0_value,
                last_k_mer_hash_strand_1_value,
                seq[k_mer_size + i],
            )
            k_mer_hashes_strand_0.add(new_k_mer_hash_strand_0_value)
            k_mer_hashes_strand_1.add(new_k_mer_hash_strand_1_value)
        else:
            (
                k_mer_hash_strand_0_index,
                k_mer_hash_strand_0_value,
            ) = k_mer_hashes_strand_0.min
            (
                k_mer_hash_strand_1_index,
                k_mer_hash_strand_1_value,
            ) = k_mer_hashes_strand_1.min
            if k_mer_hash_strand_0_value <= k_mer_hash_strand_1_value:
                minimizers_set.add((k_mer_hash_strand_0_index, 0))
                # yield k_mer_hash_strand_0_value, 0, 0
                yield k_mer_hash_strand_0_value, k_mer_hash_strand_0_index, 0
            else:
                minimizers_set.add((k_mer_hash_strand_1_index, 1))
                # yield k_mer_hash_strand_1_value, 0, 1
                yield k_mer_hash_strand_1_value, k_mer_hash_strand_1_index, 1

        for j in range(seq_length - k_mer_size - window_size + 1):
            last_k_mer_hash_strand_0_value = k_mer_hashes_strand_0.last_value
            last_k_mer_hash_strand_1_value = k_mer_hashes_strand_1.last_value
            (
                new_k_mer_hash_strand_0_value,
                new_k_mer_hash_strand_1_value,
            ) = self.calc_hash(
                last_k_mer_hash_strand_0_value,
                last_k_mer_hash_strand_1_value,
                seq[k_mer_size + window_size + j - 1],
            )
            (
                k_mer_hash_strand_0_index,
                k_mer_hash_strand_0_value,
            ) = k_mer_hashes_strand_0.add(new_k_mer_hash_strand_0_value)

            (
                k_mer_hash_strand_1_index,
                k_mer_hash_strand_1_value,
            ) = k_mer_hashes_strand_1.add(new_k_mer_hash_strand_1_value)
            if (
                k_mer_hash_strand_0_value <= k_mer_hash_strand_1_value
                and (k_mer_hash_strand_0_index, 0) not in minimizers_set
            ):
                minimizers_set.add((k_mer_hash_strand_0_index, 0))
                # yield k_mer_hash_strand_0_value, j + 1, 0
                yield k_mer_hash_strand_0_value, k_mer_hash_strand_0_index, 0
            if (
                k_mer_hash_strand_0_value > k_mer_hash_strand_1_value
                and (k_mer_hash_strand_1_index, 1) not in minimizers_set
            ):
                minimizers_set.add((k_mer_hash_strand_1_index, 1))
                # yield k_mer_hash_strand_1_value, j + 1, 1
                yield k_mer_hash_strand_1_value, k_mer_hash_strand_0_index, 1

    def run(self) -> None:
        ref_seq_dict: dict[RefId, str] = self.parse_ref_seq_file()
        read_seq_dict: dict[ReadId, str] = self.parse_read_seq_file()

        ref_minimizer_dict: collections.defaultdict[
            Minimizer, set[tuple[RefId, Index, Strand]]
        ] = collections.defaultdict(set)
        for ref_id, ref_seq in ref_seq_dict.items():
            for (minimizer, ref_pos, strand) in self.sketch_minimizers(ref_seq):
                ref_minimizer_dict[minimizer].add((ref_id, ref_pos, strand))
        # print(ref_minimizer_dict)

        outputs = []
        for read_id, read_seq in read_seq_dict.items():
            # hits: dict[tuple[RefId, Strand], Index] = {}
            # hits: collections.defaultdict[
            #     tuple[RefId, Strand], list[tuple[Index, Index]]
            # ] = collections.defaultdict(list)
            read_seq_length = len(read_seq)
            overlap_index_difference = (
                read_seq_length - self.window_size - self.k_mer_size + 1
            )
            hits: collections.defaultdict[
                RefId,
                collections.Counter,
            ] = collections.defaultdict(collections.Counter)
            for (read_minimizer, read_pos, read_strand) in self.sketch_minimizers(read_seq):
                for (ref_id, ref_pos, ref_strand) in ref_minimizer_dict[read_minimizer]:
                    # print(read_minimizer, read_pos, read_strand, ref_id, ref_pos, ref_strand)
                    if read_strand == ref_strand:
                        hits[ref_id][ref_pos - read_pos] += 1
                        # hits[(ref_id, 0)] = min(ref_pos-read_pos, hits.get((ref_id, 0), ref_pos-read_pos))
                        # hits[(ref_id, 0)].append((ref_pos, read_pos))
                    else:
                        hits[ref_id][
                            -(ref_pos + read_pos - overlap_index_difference)
                        ] += 1
                        # hits[(ref_id, 1)] = min(ref_pos-read_pos, hits.get((ref_id, 1), ref_pos-read_pos))
                        # hits[(ref_id, 1)].append((ref_pos, read_pos))
            else:
                score = 0
                output = ""
                for ref_id, counter in hits.items():
                    most_common_hit = counter.most_common(1)[0]
                    if most_common_hit[1] > score:
                        score = most_common_hit[1]
                        pos = (
                            most_common_hit[0] + 1
                            if most_common_hit[0] >= 0
                            else -most_common_hit[0] + 1
                        )
                        output = self.format(
                            read_id, ref_id, pos, int(most_common_hit[0] < 0)
                        )
                else:
                    outputs.append(output)
                # hits[(ref_id, 1)] = min(ref_pos*2-read_pos, hits.get((ref_id, 1), ref_pos-read_pos))

            # print(read_id)
            # print(hits)
            # print()
            # for (ref_id, strand_order), ref_pos in hits.items():
            #     if strand_order == 0 and read_seq[0] == ref_seq_dict[ref_id][ref_pos]:
            #         outputs.append(self.format(read_id, ref_id, ref_pos+1, 0))
            #         break
            #     if strand_order == 1 and self.seq2hash(read_seq[-1]) == self.seq2hash(ref_seq_dict[ref_id][ref_pos-1], strand=1):
            #         outputs.append(self.format(read_id, ref_id, ref_pos, 1))
            #         break
        with self.mapping_result_output_file_path.open("w") as f:
            f.writelines(outputs)

    def format(self, read_id: ReadId, ref_id: RefId, pos: int, strand_order: Strand):
        return f'{read_id}\t{ref_id}\t{pos}\t{"+-"[strand_order]}\n'

    @cache
    def calc_hash(
        self,
        last_k_mer_hash_strand_0_value: Minimizer,
        last_k_mer_hash_strand_1_value: Minimizer,
        next_nucleotide: Literal["A", "C", "G", "T"],
    ) -> tuple[Minimizer, Minimizer]:
        k_mer_hash_mod = self.k_mer_hash_mod
        k_mer_hash_strand_0_value = (
            last_k_mer_hash_strand_0_value % k_mer_hash_mod
        ) * 4 + self.seq2hash(next_nucleotide)
        k_mer_hash_strand_1_value = (
            last_k_mer_hash_strand_1_value >> 2
        ) + self.seq2hash(next_nucleotide, strand=1) * k_mer_hash_mod
        return k_mer_hash_strand_0_value, k_mer_hash_strand_1_value

    @staticmethod
    @cache
    def seq2hash(seq: str, *, strand: Strand = 0) -> int:
        res: int = 0
        if strand == 0:
            for i, v in enumerate(seq[::-1]):
                res += NUCLEOTIDE_TO_INTEGER_MAPPING[v] * 4**i
        else:
            for i, v in enumerate(seq[::-1]):
                res += (3 - NUCLEOTIDE_TO_INTEGER_MAPPING[v]) * 4**i
        return res


if __name__ == "__main__":
    IS_TEST = False
    if IS_TEST:
        minimap = Minimap(
            CURRENT_FOLDER_DIR / "test" / "read.fastq",
            CURRENT_FOLDER_DIR / "test" / "ref.fasta",
            window_size=10,
            k_mer_size=20,
        )
        minimap.run()
    else:
        minimap = Minimap(
            CURRENT_FOLDER_DIR.parent / "SE11" / "Illumina_SE11.fastq",
            CURRENT_FOLDER_DIR.parent / "SE11" / "ref_SE11.fasta",
            window_size=10,
            k_mer_size=20,
        )
        minimap.run()
