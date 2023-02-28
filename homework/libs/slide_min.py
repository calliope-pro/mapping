from collections import deque
from typing import cast, TypeAlias, TypeVar

Index: TypeAlias = int
Value = TypeVar('Value', int , str)


class SlideMin:
    deq: deque[tuple[Index, Value]]
    section_length: int

    def __init__(self, init_value: Value, section_length: int) -> None:
        self.deq = deque([(0, init_value)])
        self.section_length = section_length

    def add(self, new_value: Value) -> tuple[Index, Value]:
        new_index: Index = cast(Index, self.deq[-1][0] + 1)
        while len(self.deq) and self.deq[-1][1] >= new_value:
            self.deq.pop()
        if len(self.deq) and self.deq[0][0] <= new_index - self.section_length:
            self.deq.popleft()
        self.deq.append((new_index, new_value))
        return self.min

    @property
    def last_value(self):
        return self.deq[-1][1]

    @property
    def min(self):
        return self.deq[0]


if __name__ == "__main__":
    a = ["AC", "AA", "CG", "CT", "GC", "AA"]
    s = SlideMin("TT", 2)
    for v in a:
        print(s.deq)
        print(s.add(v))
