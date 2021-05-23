# Minimal example of using BioC classes in Webassembly

First run this to generate the WASM and JS module:

```sh
emcc --bind -I../include -std=c++1z target.cpp -o target.js
```

Create the mock files with:

```r
Rdevel -f generator.R
```

Then set up a local HTTP server:

```sh
python3 -m http.server
```

And open up `target.html`.
On the console you should see something like `"SE dimensions: 26 10"`.

