To make a dot file of the classes structure, use pyreverse:

```sh
cd cpg_pipes
pyreverse .
```

Manually customise `classes.dot`, then render with

```sh
dot -Tpng classes.dot > classes.png
```
