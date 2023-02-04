import collections
from datetime import datetime
from pathlib import Path
from time import perf_counter
from typing import Literal, NewType

Minimizer = NewType("Minimizer", int)
RefId = NewType("RefId", str)
ReadId = NewType("ReadId", str)
Strand = Literal[0, 1]

CURRENT_FOLDER_DIR = Path(__file__).resolve().parent

NUCLEOTIDE_TO_INTEGER_MAPPING: dict[
    Literal["A", "C", "G", "T"], Literal[0, 1, 2, 3]
] = {"A": 0, "C": 1, "G": 2, "T": 3}


class Minimap:
    read_seq_file_path: Path
    ref_seq_file_path: Path
    mapping_result_output_file_path: Path

    def __init__(
        self,
        read_seq_file_path: Path,
        ref_seq_file_path: Path,
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

    def sketch_minimizers(
        self, seq: str, window_size: int, k_mer_size: int
    ) -> set[tuple[Minimizer, int, Strand]]:
        assert (
            window_size <= k_mer_size
        ), "k_mer_size must be greater than or equal to window_size"
        minimizers_set: set[tuple[Minimizer, int, Strand]] = set()
        seq_length = len(seq)
        k_mer = seq[:k_mer_size]
        k_mer_hashes = collections.deque([self.seq2hash(k_mer)])
        for i in range(window_size - 1):
            k_mer_hashes.append(
                (k_mer_hashes[-1] % (4 ** (k_mer_size - 1))) * 4
                + self.seq2hash(seq[k_mer_size + i])
            )

        for i in range(seq_length - k_mer_size - window_size + 1):
            minimizer = Minimizer(min(k_mer_hashes))
            minimizers_set.add((minimizer, i, 1))
            k_mer_hashes.append(
                (k_mer_hashes[-1] % (4 ** (k_mer_size - 1))) * 4
                + self.seq2hash(seq[k_mer_size + window_size - 1 + i])
            )
            k_mer_hashes.popleft()
        else:
            minimizer = Minimizer(min(k_mer_hashes))
            minimizers_set.add(
                (minimizer, seq_length - k_mer_size - window_size + 1, 1)
            )
        return minimizers_set

    def run(self, hits_length_threshold=4) -> None:
        ref_seq_dict: dict[RefId, str] = self.parse_ref_seq_file()
        ref_seq_dict_keys = list(ref_seq_dict.keys())
        read_seq_dict: dict[ReadId, str] = self.parse_read_seq_file()

        ref_minimizer_dict: collections.defaultdict[
            Minimizer, set[tuple[int, int, Strand]]
        ] = collections.defaultdict(lambda: set())
        for ref_seq_dict_idx, ref_seq in enumerate(ref_seq_dict.values()):
            minimizers_set = self.sketch_minimizers(ref_seq, 20, 20)
            for (minimizer, pos, strand) in minimizers_set:
                ref_minimizer_dict[minimizer].add((ref_seq_dict_idx, pos, strand))

        outputs = []
        for read_id, read_seq in read_seq_dict.items():
            hits: list[tuple[int, Strand, int, int]] = []
            minimizers_set = self.sketch_minimizers(read_seq, 20, 20)
            for (read_minimizer, read_pos, read_strand) in minimizers_set:
                for (ref_seq_dict_idx, ref_pos, ref_strand) in ref_minimizer_dict[
                    read_minimizer
                ]:
                    if read_strand == ref_strand:
                        hits.append((ref_seq_dict_idx, 0, read_pos - ref_pos, ref_pos))
                    else:
                        hits.append((ref_seq_dict_idx, 1, read_pos + ref_pos, ref_pos))
            hits.sort()

            left_idx = 0
            for i in range(len(hits)):
                if (
                    i == len(hits) - 1
                    or hits[i][0] != hits[i + 1][0]
                    or hits[i][1] != hits[i + 1][1]
                    or hits[i + 1][2] - hits[i][2] >= hits_length_threshold
                ):
                    mapped = hits[left_idx : i + 1]
                    # mapped_counter = collections.Counter(map(lambda x: x[-1], mapped)).most_common()
                    # mapped_counter.sort(key=lambda x: (-x[1], x[0]))
                    outputs.append(self.format(read_id, ref_seq_dict_keys[mapped[0][0]], mapped[0][2], mapped[0][1]))
                    left_idx = i + 1
                    break
        with self.mapping_result_output_file_path.open('w') as f:
            f.writelines(outputs)

    def format(self, read_id: ReadId, ref_id: RefId, pos: int, strand: Strand):
        return f'{read_id}\t{ref_id}\t{pos}\t{"+-"[strand]}\n'

    @staticmethod
    def seq2hash(seq: str) -> int:
        res: int = 0
        for i, v in enumerate(seq[::-1]):
            res += NUCLEOTIDE_TO_INTEGER_MAPPING[v] * 4**i
        return res


if __name__ == "__main__":
    IS_TEST = False
    if IS_TEST:
        minimap = Minimap(
            CURRENT_FOLDER_DIR / "test" / "read.fastq",
            CURRENT_FOLDER_DIR / "test" / "ref.fasta",
        )
        minimap.run()
    else:
        minimap = Minimap(
            CURRENT_FOLDER_DIR.parent / "SE11" / "Illumina_SE11.fastq",
            CURRENT_FOLDER_DIR.parent / "SE11" / "ref_SE11.fasta",
        )
        minimap.run()

