This directory implements the cross-subvector aggregation algorithm of HLE-DVC.
The implementation in this directory is a modified version of Hyperproofs-Batch
(https://github.com/hyperproofs/hyperproofs
 and https://github.com/hyperproofs/gipa-go).

Below we provide a small demo. One can run the following commands from the project root directory:
```text
go run main.go
go test -v ./batch -bench=. -run=Bench -benchtime 8x -timeout 240m
```

The output is as follows:


```text
goos: linux
goarch: amd64
pkg: github.com/hyperproofs/gipa-go/batch
cpu: Intel(R) Core(TM) i5-14600KF
BenchmarkBatch
BenchmarkBatch/4/Prove;64
BenchmarkBatch/4/Prove;64-12                   2         106187184 ns/op
BenchmarkBatch/4/Verifier;64
BenchmarkBatch/4/Verifier;64-12                2          27381707 ns/op
PASS
ok      github.com/hyperproofs/gipa-go/batch    0.423s
```

To reproduce our experimental results, please set
const `MAX_AGG_SIZE = 1 << 19` in `main.go`, and set
`M := uint32(1024)` in `batch/batch_bench_test.go`.
