from typing import Optional

from pydantic import BaseModel


class ChromatinRegion(BaseModel):
    chromosome: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start

    @staticmethod
    def from_string(region: str) -> "ChromatinRegion":
        chromosome, start_end = region.split(":")
        start, end = start_end.split("-")

        if not chromosome.startswith("chr"):
            chromosome = f"chr{chromosome}"

        try:
            start = int(start)
        except ValueError:
            raise ValueError(f"Invalid start position: {start}")

        try:
            end = int(end)
        except ValueError:
            raise ValueError(f"Invalid end position: {end}")

        if start < 0:
            raise ValueError(f"Start position must be greater than 0: {start}")

        if end < 0:
            raise ValueError(f"End position must be greater than 0: {end}")

        if start >= end:
            raise ValueError(f"Start position must be less than end position: {start} >= {end}")

        return ChromatinRegion(
            chromosome=chromosome,
            start=start,
            end=end
        )

    def __str__(self) -> str:
        return f"{self.chromosome}:{self.start}-{self.end}"

    def __repr__(self) -> str:
        return f"{self.chromosome}:{self.start}-{self.end}"

    def __len__(self) -> int:
        return self.end - self.start

    def __contains__(self, other: "ChromatinRegion") -> bool:
        return self.start <= other.start and self.end >= other.end

    def __add__(self, other: "ChromatinRegion") -> "ChromatinRegion":
        if self.chromosome != other.chromosome:
            raise ValueError("Cannot add regions from different chromosomes")

        return ChromatinRegion(
            chromosome=self.chromosome,
            start=min(self.start, other.start),
            end=max(self.end, other.end)
        )

    def __sub__(self, other: "ChromatinRegion") -> Optional["ChromatinRegion"]:
        if self.chromosome != other.chromosome:
            raise ValueError("Cannot subtract regions from different chromosomes")

        if self.start >= other.start and self.end <= other.end:
            return None

        if self.start >= other.start and self.end >= other.end:
            return ChromatinRegion(
                chromosome=self.chromosome,
                start=other.end,
                end=self.end
            )

        if self.start <= other.start and self.end <= other.end:
            return ChromatinRegion(
                chromosome=self.chromosome,
                start=self.start,
                end=other.start
            )

        return ChromatinRegion(
            chromosome=self.chromosome,
            start=self.start,
            end=self.end
        )


class Enhancer3dProject(BaseModel):
    id: str
    author: Optional[str] = None
    species: str
    cell_line: str


class Enhancer3dProjectDataset(BaseModel):
    ensemble_id: str
    ensemble_region: ChromatinRegion


class Enhancer3dProjectConfiguration(BaseModel):
    enhancer_atlas_dataset_name: str
    gencode_annotation_dataset_name: str
    base_pair_linear_distance_threshold: Optional[int] = None
    enhancer_promoter_pairs_chunk_size: int = 10000
