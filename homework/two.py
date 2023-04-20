import collections
from collections.abc import Iterator
from datetime import datetime
from functools import cache
from multiprocessing import Pool, cpu_count
from multiprocessing.pool import Pool as IPool
from pathlib import Path
from typing import Literal, TypeAlias

from libs.slide_min import Index, SlideMin, Value

Minimizer: TypeAlias = int
RefId: TypeAlias = str
ReadId: TypeAlias = str
Strand = Literal[0, 1]
Twist = Literal[0, 1]

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
    processes: int

    def __init__(
        self,
        read_seq_file_path: Path,
        ref_seq_file_path: Path,
        window_size: int,
        k_mer_size: int,
        mapping_result_output_file_path: Path | None = None,
        processes: int | None = None,
    ) -> None:
        self.read_seq_file_path = read_seq_file_path
        self.ref_seq_file_path = ref_seq_file_path
        self.mapping_result_output_file_path = (
            CURRENT_FOLDER_DIR
            / f'mapping_result_PacBio_{datetime.now().isoformat(timespec="seconds")}.tsv'
            if mapping_result_output_file_path is None
            else mapping_result_output_file_path
        )
        self.k_mer_size = k_mer_size
        self.k_mer_hash_mod = (
            1 << ((k_mer_size - 1) * 2)
        ) - 1  # = 4 ** (k_mer_size - 1) - 1, 下位bit全て1
        self.window_size = window_size
        self.processes = (
            cpu_count() if processes is None else min(processes, cpu_count())
        )
        print(f"This program is runnig with {self.processes} processes.")

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
        minimizers_set: set[Index] = set()
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
            if k_mer_hash_strand_0_value < k_mer_hash_strand_1_value:
                minimizers_set.add(k_mer_hash_strand_0_index)
                yield k_mer_hash_strand_0_value, k_mer_hash_strand_0_index, 0
            if k_mer_hash_strand_0_value > k_mer_hash_strand_1_value:
                minimizers_set.add(k_mer_hash_strand_1_index)
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
                k_mer_hash_strand_0_value < k_mer_hash_strand_1_value
                and k_mer_hash_strand_0_index not in minimizers_set
            ):
                minimizers_set.add(k_mer_hash_strand_0_index)
                yield k_mer_hash_strand_0_value, k_mer_hash_strand_0_index, 0
            if (
                k_mer_hash_strand_0_value > k_mer_hash_strand_1_value
                and k_mer_hash_strand_1_index not in minimizers_set
            ):
                minimizers_set.add(k_mer_hash_strand_1_index)
                yield k_mer_hash_strand_1_value, k_mer_hash_strand_1_index, 1

    def generate_analyzed_read_seq_output(
        self,
        args: tuple[
            ReadId,
            str,
            collections.defaultdict[Minimizer, list[tuple[RefId, Index, Strand]]],
        ],
    ) -> str:
        read_id, read_seq, ref_minimizer_dict = args
        read_seq_length = len(read_seq)
        hits: collections.defaultdict[
            RefId, collections.Counter
        ] = collections.defaultdict(collections.Counter)
        read_minimizer_iterator = self.sketch_minimizers(read_seq)
        for (read_minimizer, read_pos, read_strand) in read_minimizer_iterator:
            for (ref_id, ref_pos, ref_strand) in ref_minimizer_dict[read_minimizer]:
                if read_strand == ref_strand:
                    hits[ref_id][ref_pos - read_pos + 1] += 1
                else:
                    hits[ref_id][
                        -(ref_pos - (read_seq_length - self.k_mer_size - read_pos) + 1)
                    ] += 1

        max_cnt = 0
        output = ""
        for ref_id, counter in hits.items():
            most_common_pos, most_common_cnt = counter.most_common(1)[0]
            if most_common_cnt > max_cnt:
                max_cnt = most_common_cnt
                output = self.format(
                    read_id, ref_id, most_common_pos, int(most_common_pos < 0)
                )
        return output

    def run(self) -> None:
        ref_seq_dict: dict[RefId, str] = self.parse_ref_seq_file()
        read_seq_dict: dict[ReadId, str] = self.parse_read_seq_file()

        ref_minimizer_dict: collections.defaultdict[
            Minimizer, list[tuple[RefId, Index, Strand]]
        ] = collections.defaultdict(list)
        for ref_id, ref_seq in ref_seq_dict.items():
            ref_minimizer_iterator = self.sketch_minimizers(ref_seq)
            for (minimizer, ref_pos, strand) in ref_minimizer_iterator:
                ref_minimizer_dict[minimizer].append((ref_id, ref_pos, strand))

        outputs: list[str] = []
        with Pool(processes=self.processes) as pool:
            pool: IPool
            outputs: list[str] = pool.map(
                self.generate_analyzed_read_seq_output,
                (
                    (read_id, read_seq, ref_minimizer_dict)
                    for read_id, read_seq in read_seq_dict.items()
                ),
            )

        with self.mapping_result_output_file_path.open("w") as f:
            f.writelines(outputs)

    def format(self, read_id: ReadId, ref_id: RefId, pos: int, strand_order: Twist):
        return f'{read_id}\t{ref_id}\t{abs(pos)}\t{"+-"[strand_order]}\n'

    def calc_hash(
        self,
        last_k_mer_hash_strand_0_value: Minimizer,
        last_k_mer_hash_strand_1_value: Minimizer,
        next_nucleotide: Literal["A", "C", "G", "T"],
    ) -> tuple[Minimizer, Minimizer]:
        k_mer_hash_mod = self.k_mer_hash_mod
        k_mer_hash_strand_0_value = (
            (last_k_mer_hash_strand_0_value & k_mer_hash_mod) << 2
        ) | self.seq2hash(next_nucleotide)
        k_mer_hash_strand_1_value = (last_k_mer_hash_strand_1_value >> 2) | (
            self.seq2hash(next_nucleotide, strand=1) << (2 * (self.k_mer_size - 1))
        )
        return k_mer_hash_strand_0_value, k_mer_hash_strand_1_value

    @staticmethod
    @cache
    def seq2hash(seq: str, *, strand: Strand = 0) -> int:
        res: int = 0
        if strand == 0:
            for v in seq:
                res = (res << 2) | NUCLEOTIDE_TO_INTEGER_MAPPING[v]
        else:
            for v in seq[::-1]:
                res = (res << 2) | (3 - NUCLEOTIDE_TO_INTEGER_MAPPING[v])
        return res


if __name__ == "__main__":
    IS_TEST = False
    if IS_TEST:
        minimap = Minimap(
            CURRENT_FOLDER_DIR / "test" / "PacBio.fastq",
            CURRENT_FOLDER_DIR / "test" / "ref.fasta",
            window_size=10,
            k_mer_size=20,
        )
        minimap.run()

    else:
        minimap = Minimap(
            CURRENT_FOLDER_DIR.parent / "SE11" / "PacBio_SE11.fastq",
            CURRENT_FOLDER_DIR.parent / "SE11" / "ref_SE11.fasta",
            window_size=10,
            k_mer_size=20,
        )
        minimap.run()
