import argparse
from enum import Enum
from pathlib import Path
from typing import List, Tuple
import spkm
import pandas as pd
from pandas import DataFrame
import numpy as np

MAX_ITER = 300
EPSILON = 0.0


class Goal(Enum):
    WEIGHT_MATRIX = "wam"
    DIAGONAL_DEGREE_MATRIX = "ddg"
    NORMALIZED_GRAPH_LAPLACIAN = "lnorm"
    JACOBI = "jacobi"
    SPKMEANS = "spk"

    def __str__(self) -> str:
        return self.value


class InvalidInputParser(argparse.ArgumentParser):
    def error(self, message: str):
        if message.startswith("unrecognized arguments:"):
            return
        print("Invalid Input!")
        exit(1)


def valid_file(path: str) -> Path:
    parsed = Path(path)
    if parsed.suffix not in [".csv", ".txt"]:
        raise TypeError("Invalid Input!")
    return parsed


parser = InvalidInputParser(description="Implementation of spkmeans algorithm")
parser.add_argument(
    "k",
    help="Number of required cluster. if equal 0 - use the eigengap heuristic",
    type=int,
)
parser.add_argument(
    "goal",
    help="Goal of the algorithm.",
    choices=list(Goal),
    type=Goal,
)
parser.add_argument(
    "file_name", help="Input file name, must end with .txt or .csv", type=valid_file
)

args = parser.parse_args()

np.random.seed(0)


def read_matrix(matrix_path: Path) -> List[List[float]]:
    with matrix_path.open("r") as f:
        lines = f.readlines()
        return [[float(num) for num in line.split(",")] for line in lines]


def read_points(points_path: Path) -> List[Tuple]:
    with points_path.open("r") as f:
        lines = f.readlines()
        return [tuple([float(val) for val in line.split(",")]) for line in lines]


def format_matrix(matrix: List[List[float]]) -> str:
    return "\n".join([format_values(row) for row in matrix])


def format_values(values: List[float]) -> str:
    return ",".join([f"{val:.4f}" for val in values])


def kmeanspp(df: DataFrame, k: int) -> List[int]:
    """Running KMeans++ algorithm"""
    vectors = df.to_numpy()
    cent_ix = [int(np.random.choice(df.index))]
    df["distance"] = 0
    for i in range(1, k):
        centroids = df.loc[cent_ix]
        for j, vec in enumerate(vectors):
            min_d = (
                np.linalg.norm(vec - centroids.drop(["distance"], axis=1), axis=1) ** 2
            )
            df.at[j, "distance"] = min(min_d)
        total_distance = sum(df["distance"])
        probs = [distance / total_distance for distance in df["distance"].tolist()]
        selected_centroid = np.random.choice(len(df), 1, p=probs)
        cent_ix.append(selected_centroid[0])
    return cent_ix


if __name__ == "__main__":
    if args.goal != Goal.JACOBI:
        points = read_points(args.file_name)
        if args.goal == Goal.WEIGHT_MATRIX:
            result = spkm.wam(points, args.k)
        elif args.goal == Goal.DIAGONAL_DEGREE_MATRIX:
            result = spkm.ddg(points, args.k)
        elif args.goal == Goal.NORMALIZED_GRAPH_LAPLACIAN:
            result = spkm.lnorm(points, args.k)
        elif args.goal == Goal.SPKMEANS:
            result = spkm.spk(points, args.k)
            k = args.k if args.k != 0 else len(result[0])
            df = pd.DataFrame(result)
            res = kmeanspp(df, k)
            print(",".join([str(i) for i in res]))
            centroids = [tuple(row[:-1]) for row in df.loc[res].to_numpy()]
            data_points = [tuple(row[:-1]) for row in df.to_numpy()]
            result = spkm.kmeans_fit(centroids, data_points, MAX_ITER, EPSILON)
        print(format_matrix(result))
    else:
        mat = read_matrix(args.file_name)
        result = spkm.jacobi(mat)
        print(format_values(result[0]))
        print(format_matrix(result[1]))
