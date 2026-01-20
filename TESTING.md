# Testing
If you have uv installed, the easiest way to test your code is to run the following from the root directory of PlantClusterFinder.py

```shell
uv sync
uv run pytest
```

## Installing rye
The [uv documentation](https://docs.astral.sh/uv/#installation) has a guide on installing uv,
but in Linux, running the following command in your home directory should install uv and set up shims.

```bash
curl -sSf https://astral.sh/uv/install.sh | bash
```
