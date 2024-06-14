#!/usr/bin/env python3

import os
import json
import sys
from typing import List, Callable
from pylatex import Document, Command, TikZ, Axis, Plot, NoEscape


def percentage_change(x: float, y: float):
    return x / y


def percentage_usage(x: float, y: float):
    return (y * 100) / (x + y)


def groupby(func: Callable[[object], object], data: List[object]):
    return (
        (key, [x for x in data if func(x) == key])
        for key in sorted(set(map(func, data)))
    )


def create_figure(filename: str, options: List[str], *plots: Plot):
    doc = Document(documentclass="standalone", document_options="tikz,border=5pt")
    doc.append(Command("usepgfplotslibrary", "colorbrewer"))
    doc.append(Command("pgfplotsset", NoEscape("colormap/Paired-9")))
    with doc.create(TikZ()) as fig:
        fig_options = [
            "width=9cm",
            "height=5.5cm",
            "cycle list name=Paired-9",
            "legend style={font=\\footnotesize}",
            *options,
        ]
        with fig.create(Axis(NoEscape(",".join(fig_options)))) as axis:
            for plot in plots:
                axis.append(plot)
    doc.generate_tex(filename)
    os.system(f"latex {filename}.tex && dvips {filename}.dvi -o {filename}.eps")


with open(sys.argv[1]) as fp:
    FULL_DATA = json.load(fp)

for (instance, n), data in groupby(
    lambda x: (x["params"]["instance"], x["params"]["n"]),
    [x for x in FULL_DATA if x["params"]["exp"] == "eval"],
):
    for i, (F, data) in enumerate(groupby(lambda x: x["params"]["F"], data)):
        create_figure(
            f"eval-{instance}-{F}".replace(".", "-"),
            [
                "xlabel={$r$ (\\% of $n$)}",
                "ylabel={speedup}",
                "ymode={log}",
                "xmin=1",
                "xmax=99",
                # "ymin=-120",
                # "ymax=320",
                # "ytick={-100, 0, 100, 200, 300}",
                "legend cell align={left}",
                "legend columns=1",
                "legend pos={north east}",
            ],
            *[
                Plot(
                    name=Command("texttt", strategy_name),
                    coordinates=[
                        (
                            (r * 100 // n),
                            percentage_change(
                                next(
                                    x["result"]["avg"]
                                    for x in entries
                                    if x["params"]["eval"] == 0
                                ),
                                next(
                                    x["result"]["avg"]
                                    for x in entries
                                    if x["params"]["eval"] == strategy_number
                                ),
                            ),
                        )
                        for r, entries in groupby(lambda x: x["params"]["r"], data)
                    ],
                )
                for strategy_number, strategy_name in list(
                    enumerate(["B", "F-INC", "R", "RV"])
                )[1:]
            ],
            Plot(
                options="mark=none,dashed,thick,color=red",
                coordinates=[(0, 1), (100, 1)],
            ),
        )

for (instance, n), data in groupby(
    lambda x: (x["params"]["instance"], x["params"]["n"]),
    [x for x in FULL_DATA if x["params"]["exp"] == "ls"],
):
    for i, (F, data) in enumerate(groupby(lambda x: x["params"]["F"], data)):
        create_figure(
            f"ls-{instance}-{F}".replace(".", "-"),
            [
                "xlabel={$r$ (\\% of $n$)}",
                "ylabel={speedup}",
                "ymode={log}",
                "xmin=1",
                "xmax=99",
                # "ymin=-120",
                # "ymax=320",
                "ytick={1,10,100,1000}",
                "legend cell align={left}",
                "legend columns=1",
                "legend pos={north east}",
            ],
            *[
                Plot(
                    name=Command("texttt", strategy_name),
                    coordinates=[
                        (
                            (r * 100 // n),
                            percentage_change(
                                next(
                                    x["result"]["dt"]
                                    for x in entries
                                    if x["params"]["eval"] == 0
                                ),
                                next(
                                    x["result"]["dt"]
                                    for x in entries
                                    if x["params"]["eval"] == strategy_number
                                ),
                            ),
                        )
                        for r, entries in groupby(lambda x: x["params"]["r"], data)
                    ],
                )
                for strategy_number, strategy_name in list(
                    enumerate(["B", "F-INC", "R", "RV"])
                )[1:]
            ],
            Plot(
                options="mark=none,dashed,thick,color=red",
                coordinates=[(0, 1), (100, 1)],
            ),
        )


os.system("rm -f *.aux *.dvi *.log")
