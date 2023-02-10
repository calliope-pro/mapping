from collections import deque
from typing import NewType

Index = NewType("Index", int)
Value = NewType("Value", int)


class SlideMin:
    deq: deque[tuple[Index, Value]]
    section_length: int

    def __init__(self, init_value: int, section_length: int) -> None:
        self.deq = deque([(Index(0), Value(init_value))])
        self.section_length = section_length

    def add(self, new_value: int) -> tuple[Index, Value]:
        new_index = self.deq[-1][0] + 1
        while len(self.deq) and self.deq[-1][1] >= new_value:
            self.deq.pop()
        if len(self.deq) and self.deq[0][0] <= new_index - self.section_length:
            self.deq.popleft()
        self.deq.append((Index(new_index), Value(new_value)))
        return self.deq[0]


if __name__ == "__main__":
    a = ['AC', 'AA', 'CG', 'CT', 'GC', 'AA']
    s = SlideMin(33, 1)
    for v in a:
        print(s.add(v))
