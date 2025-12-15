\ 実行不能問題 - IIS (Irreducible Infeasible Subsystem) のデモ用
\
\ この問題は矛盾する制約を含んでおり、解が存在しない
\ SCIPはこれを検出し、IISを計算できる

Minimize
 obj: x + y + z

Subject To
 \ 制約1: x + y >= 10
 c1: x + y >= 10

 \ 制約2: x + y <= 5
 c2: x + y <= 5

 \ 制約3: z >= 0 (これはIISに含まれない)
 c3: z >= 0

 \ 制約4: x >= 3
 c4: x >= 3

 \ 制約5: y >= 2
 c5: y >= 2

Bounds
 0 <= x <= 100
 0 <= y <= 100
 0 <= z <= 100

End
