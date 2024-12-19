from typing import List, Any

import pandas as pd


class DataFrameBufferedSink:
    __slots__ = ['columns', 'batch_size', 'index', 'buffer', 'dfs']

    def __init__(self, columns, batch_size=1000):
        self.columns: List[str] = columns
        self.batch_size: int = batch_size

        self.index: int = 0
        self.buffer: List[List[Any]] = []

        self.dfs: List[pd.DataFrame] = []

    @property
    def df(self):
        self.flush()
        return pd.concat(self.dfs, ignore_index=True)

    def flush(self):
        if self.buffer:
            # print(f'Flushing {len(self.buffer)} rows, current df count: {len(self.dfs)}')
            self.dfs.append(pd.DataFrame(self.buffer, columns=self.columns))

            del self.buffer
            self.buffer = []
            self.index = 0

    def write(self, data):
        self.index += 1
        self.buffer.append(data)
        if self.index >= self.batch_size:
            self.flush()

    def __add__(self, data):
        self.write(data)
        return self

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        del self.buffer
        del self.dfs

    def __call__(self, data):
        self.write(data)
