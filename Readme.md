# Faster Nearest Neighbor Queries on Geographic Data

This is the Code for the experiments to the Paper "Faster Nearest Neighbor Queries on Geographic Data".
Before you start, you have to get the positions of the data points as well as positions of the grid.
Running the script `get_data.sh` in the `data` directory will get the latest amenity data for the Berlin area.
We provide a nix-shell file with the appropriate python dependencies.

We provide an example configuration that can be used to run the benchmark on the Berlin data set with $r = 2.0$.
To run the example you first have to install the tool chain for the Rust programming language.
After this you need to install the wrapper for the benchmark harness using
```shell
cargo install cargo-criterion
```
Finally you can use the following command to run the example.
The first time any given configuration is run, the approximations are built and cached.
This might take a few additional minutes.
```shell
CONFIG_FILE=example.json cargo criterion --bench knn_comparison
```

You can set the benchmark time using the environment variable `RUNNING_TIME`.

## Machine Readable output
The command above gives human-readable output.
To get output that is better machine-readable use the following command.
This will create a JSON file, representing the gathered benchmark data.
```shell
export CONFIG_FILE=example.json
cargo criterion --output-format verbose --bench knn_comparison --message-format=json | tee $CONFIG_FILE.log
awk -f postprocess.awk -i $CONFIG_FILE.log
```
