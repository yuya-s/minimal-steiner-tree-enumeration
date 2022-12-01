# Cost-constrained Minimal Steiner Tree Enumeration

This is an implementation of cost-constrained steiner tree enumeration in the following paper:\
Yuya Sasaki, [Cost-constrained Minimal Steiner Tree Enumeration](https://dl.acm.org/doi/abs/10.1145/3511808.3557570)

It aims to approximately and practically enumerate minimal Steiner trees with the top-k smallest costs in a large graph by using binary decision diagram.

## setup

CMake 3.7 or later

google [densehash](https://github.com/sparsehash/sparsehash)


## dataset
Edge file includes a set of edges, terminals, and the optimal cost.

The top line has four values: the number of vertices, the number of edges, the number of terminals, and the optimal cost.
Then, it lists edges with vertex id1, vertex id 2, and its edge weight.
Finally, it lists the vertex ids as terminals.

```
1244 1971 34 1049
0 1 5
1 2 5
2 3 5
3 4 5
4 5 5
...
1238 1239 5
1239 1240 5
15
186
...
1154
1178
```

Also, we need to prepare Steiner trees as seed trees.
These files should include a set of edge ids that construct Steiner tree. We here note that edge ids indicate the line number -2 in the graph file (e.g., edge id of "0 1 5" is zero)

For example,
```
1049 // cost of Steiner tree
39
97
162
...
1892
1893
```

## run code

1. cmake CMakeLists.txt
2. make
3. ./run.sh

## arguments
Our code has the following parameters.

- i `X`: set X as the input graph
- o `X`: set X as output file name, default ./result/test
- r: graph reduction
- m `X`: the maximum BDD node size  
- k `X`: set X as the number of output trees, default 1
- th `X`: sex X as the threshold of costs
- k `X`: set X as the maximum length of paths, default 2
- s `X` `filename1` ... `filename X`: set X as the number of seed trees and filenames of seed trees

See run.sh in details.

## citing
If you find our code is useful for your research, please consider citing the following paper:

    @inproceedings{sasaki2022cost,
    title={Cost-constrained Minimal Steiner Tree Enumeration},
    author={Sasaki, Yuya},
    booktitle={Proceedings of the 31st ACM International Conference on Information \& Knowledge Management},
    pages={4439--4443},
    year={2022}
    }

## contact
Please let me know if they have problem to sasaki@ist.osaka-u.ac.jp
