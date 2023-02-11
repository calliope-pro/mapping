import collections
from collections.abc import Iterator
from datetime import datetime
from pathlib import Path
from typing import ClassVar, Literal, NewType, Type

from libs.slide_min import Index, SlideMin, Value

Minimizer = NewType("Minimizer", int)
RefId = NewType("RefId", str)
ReadId = NewType("ReadId", str)
Strand = Literal[0, 1]

CURRENT_FOLDER_DIR = Path(__file__).resolve().parent

NUCLEOTIDE_TO_INTEGER_MAPPING: dict[
    Literal["A", "C", "G", "T"], Literal[0, 1, 2, 3]
] = {"A": 0, "C": 1, "G": 2, "T": 3}


class Minimap:
    cached_seq_dict: ClassVar[
        dict[tuple[Minimizer, Literal["A", "C", "G", "T"]], tuple[Minimizer, Minimizer]]
    ]
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
        self.cached_seq_dict = (
            {} if self.cached_seq_dict is None else self.cached_seq_dict
        )
        self.k_mer_size = k_mer_size
        self.k_mer_hash_mod = 1 << ((k_mer_size - 1) * 2)
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
                    tmp_ref_id = RefId(line.split()[0][1:])
                else:
                    seq += line
            ref_seq_dict[tmp_ref_id] = seq
        return ref_seq_dict

    def parse_read_seq_file(self) -> dict[ReadId, str]:
        read_seq_dict: dict[ReadId, str] = {}
        with self.read_seq_file_path.open("r") as f:
            while True:
                read_id = ReadId(f.readline().rstrip()[1:])
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
        minimizers_set: set[tuple[Minimizer, Index, Strand]] = set()
        seq_length = len(seq)
        k_mer_init = seq[:k_mer_size]
        k_mer_hashes_strand_0 = SlideMin(self.seq2hash(k_mer_init), window_size)
        k_mer_hashes_strand_1 = SlideMin(
            self.seq2hash(k_mer_init, strand=1), window_size
        )
        for i in range(window_size - 1):
            last_k_mer_hash_strand_0_value = k_mer_hashes_strand_0.deq[-1][1]
            last_k_mer_hash_strand_1_value = k_mer_hashes_strand_1.deq[-1][1]
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
            ) = k_mer_hashes_strand_0.deq[0]
            (
                k_mer_hash_strand_1_index,
                k_mer_hash_strand_1_value,
            ) = k_mer_hashes_strand_1.deq[0]
            if k_mer_hash_strand_0_value < k_mer_hash_strand_1_value:
                minimizers_set.add(
                    (k_mer_hash_strand_0_value, k_mer_hash_strand_0_index, 0)
                )
                yield k_mer_hash_strand_0_value, k_mer_hash_strand_0_index, 0
            if k_mer_hash_strand_0_value > k_mer_hash_strand_1_value:
                minimizers_set.add(
                    (k_mer_hash_strand_1_value, k_mer_hash_strand_1_index, 1)
                )
                yield k_mer_hash_strand_1_value, k_mer_hash_strand_1_index, 1

        for j in range(seq_length - k_mer_size - window_size + 1):
            last_k_mer_hash_strand_0_value = k_mer_hashes_strand_0.deq[-1][1]
            last_k_mer_hash_strand_1_value = k_mer_hashes_strand_1.deq[-1][1]
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
                and (k_mer_hash_strand_0_value, k_mer_hash_strand_0_index, 0)
                not in minimizers_set
            ):
                minimizers_set.add(
                    (k_mer_hash_strand_0_value, k_mer_hash_strand_0_index, 0)
                )
                yield k_mer_hash_strand_0_value, k_mer_hash_strand_0_index, 0
            if (
                k_mer_hash_strand_0_value > k_mer_hash_strand_1_value
                and (k_mer_hash_strand_1_value, k_mer_hash_strand_1_index, 1)
                not in minimizers_set
            ):
                minimizers_set.add(
                    (k_mer_hash_strand_1_value, k_mer_hash_strand_1_index, 1)
                )
                yield k_mer_hash_strand_1_value, k_mer_hash_strand_1_index, 1

    def run(self, hits_length_threshold=5) -> None:
        ref_seq_dict: dict[RefId, str] = self.parse_ref_seq_file()
        read_seq_dict: dict[ReadId, str] = self.parse_read_seq_file()

        ref_minimizer_dict: collections.defaultdict[
            Minimizer, set[tuple[RefId, Index, Strand]]
        ] = collections.defaultdict(set)
        for ref_id, ref_seq in ref_seq_dict.items():
            minimizers_list = self.sketch_minimizers(ref_seq)
            for (minimizer, ref_pos, strand) in minimizers_list:
                ref_minimizer_dict[minimizer].add((ref_id, ref_pos, strand))

        outputs = []
        for read_id, read_seq in read_seq_dict.items():
            # hits: dict[tuple[RefId, Strand], Index] = {}
            hits: collections.defaultdict[
                tuple[RefId, Strand], list[tuple[Index, Index]]
            ] = collections.defaultdict(list)
            minimizers_list = self.sketch_minimizers(read_seq)
            for (read_minimizer, read_pos, read_strand) in minimizers_list:
                print(read_minimizer, read_pos, read_strand)
                for (ref_id, ref_pos, ref_strand) in ref_minimizer_dict[read_minimizer]:
                    if read_strand == ref_strand:
                        # hits[(ref_id, 0)] = min(ref_pos-read_pos, hits.get((ref_id, 0), ref_pos-read_pos))
                        hits[(ref_id, 0)].append((ref_pos, read_pos))
                    else:
                        hits[(ref_id, 1)].append((ref_pos, read_pos))
                        # hits[(ref_id, 1)] = min(ref_pos*2-read_pos, hits.get((ref_id, 1), ref_pos-read_pos))

            print(read_id)
            print(hits)
            print()
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

    def calc_hash(
        self,
        last_k_mer_hash_strand_0_value: Minimizer,
        last_k_mer_hash_strand_1_value: Minimizer,
        next_nucleotide: Literal["A", "C", "G", "T"],
    ) -> tuple[Minimizer, Minimizer]:
        k_mer_hash_mod = self.k_mer_hash_mod
        if (last_k_mer_hash_strand_0_value, next_nucleotide) in self.cached_seq_dict:
            return self.cached_seq_dict[
                (last_k_mer_hash_strand_0_value, next_nucleotide)
            ]

        k_mer_hash_strand_0_value = (
            last_k_mer_hash_strand_0_value % k_mer_hash_mod
        ) * 4 + self.seq2hash(next_nucleotide)
        k_mer_hash_strand_1_value = (
            last_k_mer_hash_strand_1_value >> 2
        ) + self.seq2hash(next_nucleotide, strand=1) * k_mer_hash_mod
        self.update_cached_seq_dict(
            (last_k_mer_hash_strand_0_value, next_nucleotide),
            (k_mer_hash_strand_0_value, k_mer_hash_strand_1_value),
        )
        return k_mer_hash_strand_0_value, k_mer_hash_strand_1_value

    def update_cached_seq_dict(
        self,
        key: tuple[Minimizer, Literal["A", "C", "G", "T"]],
        value: tuple[Minimizer, Minimizer],
    ) -> None:
        self.cached_seq_dict[key] = value

    @staticmethod
    def seq2hash(seq: str, *, strand: Strand = 0) -> int:
        res: int = 0
        if strand == 0:
            for i, v in enumerate(seq[::-1]):
                res += NUCLEOTIDE_TO_INTEGER_MAPPING[v] * 4**i
        else:
            for i, v in enumerate(seq):
                res += (3 - NUCLEOTIDE_TO_INTEGER_MAPPING[v]) * 4**i
        return res

    @classmethod
    def minimap_cls_factory(
        cls: Type["Minimap"],
        cached_seq_dict: dict[
            tuple[Minimizer, Literal["A", "C", "G", "T"]], tuple[Minimizer, Minimizer]
        ],
    ) -> Type["Minimap"]:
        cls.cached_seq_dict = cached_seq_dict
        return cls


if __name__ == "__main__":
    IS_TEST = True
    if IS_TEST:
        minimap_factory = Minimap.minimap_cls_factory({})
        minimap = minimap_factory(
            CURRENT_FOLDER_DIR / "test" / "read.fastq",
            CURRENT_FOLDER_DIR / "test" / "ref.fasta",
            window_size=1,
            k_mer_size=3,
        )
        minimap.run()
    else:
        minimap_factory = Minimap.minimap_cls_factory({})
        minimap = minimap_factory(
            CURRENT_FOLDER_DIR / "test" / "read.fastq",
            CURRENT_FOLDER_DIR / "test" / "ref.fasta",
        )
        minimap.run(window_size=4, k_mer_size=25)
