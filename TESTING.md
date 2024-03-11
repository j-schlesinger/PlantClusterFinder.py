# Testing
If you have rye installed, the easiest way to test your code is to run the following from the root directory of PlantClusterFinder.py

```shell
rye sync
rye run pytest
```

## Installing rye
The [rye documentation](https://rye-up.com/guide/installation/#installing-rye) has a guide on installing rye,
but in Linux, running the following command in your home directory should install rye and set up shims.

```bash
curl -sSf https://rye-up.com/get | bash
echo 'source "$HOME/.rye/env"' >> ~/.profile
```
