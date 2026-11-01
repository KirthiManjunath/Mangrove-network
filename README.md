# Mangrove Maze

A Complex Networks project (IIT Madras, Sem 5): **Sundarbans tidal channels as networks**.
We extract a channel network from each Sentinel-2 NDWI image of the Sundarbans (October–November, 2016–2026) using [RivGraph](https://github.com/VeinsOfTheEarth/RivGraph). Then we study how its structure compares with random models and how it changes over the years.

The project lives in two places that hold **the same folder**:

| Where | What it is for |
|---|---|
| Google Drive: `MyDrive/Mangrove_Maze` | Where we work. Every Colab notebook runs from this folder. |
| GitHub: [KirthiManjunath/Mangrove-network](https://github.com/KirthiManjunath/Mangrove-network) | Version history and backup. The git repo *is* the Drive folder (`Mangrove_Maze/.git`). |

> ⚠️ The GitHub repo is **public** and the Drive folder is link-shared. Never put tokens, passwords or personal data in a notebook, a cell output or a file in this folder.

---

## What is in the folder

| Path | What it is |
|---|---|
| `Mangroves.ipynb` | **Setup notebook.** Mounts Drive, installs packages from the wheelhouse, checks versions, and holds the git cells (connect, `.gitignore`, commit and push). It also has the early exploration: two ways to turn an image into a network, and a first look at degree distribution, betweenness and correlations. |
| `main.ipynb` | **The graphs everyone works on (unweighted).** It builds the yearly networks **G1 … G10**, a random-tree baseline **B1 … B10** for each year, saves them to `Results/timeseries/`, and plots the parameters over the years. See [The shared graphs](#the-shared-graphs-mainipynb). |
| `main_weighted.ipynb` | **The same graphs, simplified and weighted.** Each link gets `conductance` = width / length. Redundant links are contracted. Weighted parameters (resistance, efficiency, λ₂, betweenness Gini in L/w) are compared with baseline trees that have the same channel geometry. Saves to `Results/timeseries_weighted/`. See [The weighted graphs](#the-weighted-graphs-main_weightedipynb). |
| `Data/` | Input images: Sentinel-2 L2A NDWI renders, one per year, 2016–2026. Named `YYYY-MM-DD-00_00_YYYY-MM-DD-23_59_Sentinel-2_L2A_NDWI.tiff`. 2018 is in the folder but `main.ipynb` skips it. |
| `Results/` | Everything the notebooks write. `main.ipynb` writes to `Results/timeseries/`, `main_weighted.ipynb` to `Results/timeseries_weighted/`, and each person to `Results/<your-name>/`. |
| `team/TEMPLATE.ipynb` | **Starting notebook for your own work.** Copy it; don't edit it. See [Using the template](#using-the-template). |
| `team/<your-name>/` | Your own notebooks. |
| `README.md` | This file. |
| `sundarbans_rg.py` | Helper functions: image → water mask → shoreline → pruned RivGraph network. Used by `main.ipynb`. |
| `tree_nulls.py` | Random-tree baselines and the comparison metrics, unweighted and weighted. Used by both main notebooks. |
| `sundarbans_simplified_network.py` | Adds weights and contracts redundant links: `build_simplified_network`, `to_unweighted`, `load_graphs`, `assign_baseline_weights`, `count_cycles`. Used by `main_weighted.ipynb`. **Must be in the folder before `main_weighted.ipynb` can run.** |
| `mask_tools.py` | Helpers for building a water mask from an RGB or NDWI image. |
| `RivGraph/` | A copy of the RivGraph source code, installed from the wheelhouse. Do not edit it unless the whole team agrees. |
| `requirements-py313.txt` | Pinned package versions for Colab (Python 3.13). |
| `wheels_py313/` | Pre-built packages (about 185 MB), so installing takes about 30 s with no downloads. Rebuilt by `Mangroves.ipynb` step 3. It is listed in `.gitignore`, so it exists on Drive only. |

---

## First-time setup (once per person)

1. **Get access.** Ask Kirthi to share the Drive folder `Mangrove_Maze` with your Google account as an **Editor** (Viewers can't save notebooks or results). Also ask for collaborator access to the GitHub repo.
2. **Add the folder to your My Drive.** `Mangrove_Maze` is a folder in Kirthi's Drive shared with you, not a Shared Drive. When you mount Drive in Colab you see only *your* My Drive, so the folder isn't there until you add a shortcut:
   - In drive.google.com, open *Shared with me*. Right-click `Mangrove_Maze` and choose *Organize → Add shortcut*. Pick **My Drive** itself (not a subfolder) and keep the name `Mangrove_Maze`.
   - Colab follows the shortcut, so the project appears at `/content/drive/MyDrive/Mangrove_Maze`. Anything you write there goes to the real shared folder, and teammates see it.
   - If Colab can't find it: check that the shortcut wasn't renamed, or do *Runtime → Disconnect and delete runtime* and mount again (a new shortcut can take a minute to appear). The template's SETUP cell also searches one and two folders deep in case you put the shortcut in a subfolder.
   - Never use *Make a copy* on the folder. That makes a private copy nobody else sees.
3. **Make a GitHub token.**
   - Go to GitHub → *Settings → Developer settings → Personal access tokens → Tokens (classic)*. Create a token with the `public_repo` scope.
   - Use a classic token: a fine-grained token can't reach a repo owned by another user.
   - In Colab, open the key icon (Secrets) and save the token as **`GITHUB_TOKEN`**. Turn on notebook access for it. Never paste the token into a cell.
4. **Run the setup once.** Open `Mangroves.ipynb` and run **Env setup** steps 1, 4 and 5.
   - Skip steps 2 and 3. They rebuild the wheelhouse, which is only needed when versions change.
   - If step 4 prints `Uninstalling numpy/scipy/pandas`, do *Runtime → Restart session*, then run step 5.

---

## How we work

### 1. Each person works in their own notebook
- Make your notebook inside the Drive folder at **`team/<your-name>/<your-name>_<topic>.ipynb`**. For example: `team/asha/asha_robustness.ipynb`. Start from `team/TEMPLATE.ipynb` (next section).
- Write your outputs to **`Results/<your-name>/`**, so nobody overwrites anyone else's files.
- **Don't edit `Mangroves.ipynb`, `main.ipynb`, `main_weighted.ipynb` or `team/TEMPLATE.ipynb`** to do your own analysis. Change them only when the setup or the shared graphs need to change, and tell the team first.

### 2. Work on the shared graphs
Everyone analyses **the same graphs, the ones `main.ipynb` and `main_weighted.ipynb` build**:
- **Unweighted questions** (cycles, bridges, paths counted in links): use the graphs from `main.ipynb`.
- **Questions that need channel size** (flow, resistance, weighted centrality): use the graphs from `main_weighted.ipynb`.

Load the saved graphs; don't extract your own network from the images. That way every result in the project refers to the same networks.

If you think the graphs need to change (for example, different pruning or another year), raise it with the team. If the team agrees, the change goes into the main notebooks, and everyone reruns them.

---

## Using the template

`team/TEMPLATE.ipynb` already has everything a new notebook needs.

**Make your copy:**
1. Open `Mangrove_Maze/team/TEMPLATE.ipynb` in Colab.
2. *File → Save a copy in Drive*. Colab saves the copy to *your* `My Drive/Colab Notebooks/` and opens it.
3. In Drive, move the copy into `Mangrove_Maze/team/<your-name>/` (create the folder first), and rename it `<your-name>_<topic>.ipynb`. Until you move it, the copy is private and isn't in the project.
4. Open the moved copy from Drive and do all your work there.

**Fill it in:**
- Cell 2, `ME = 'your-name'`: your name. It sets your output folder `Results/<your-name>/`. The cell stops with an error until you change it.
- Last cell, `GH_USER, GH_EMAIL`: your GitHub username and email.
- Cell 6 (EXAMPLE): replace it with your own analysis. Keep cells 1–5.

**What each cell does:**

| Cell | What it does | Run it |
|---|---|---|
| 1. SETUP | Mounts Drive, finds `Mangrove_Maze` (also if the shortcut is in a subfolder), moves into it, installs packages from the wheelhouse. | Every session, first. If it prints `Uninstalling numpy/scipy/pandas`, do *Runtime → Restart session* and run it again. |
| 2. ME | Sets your name and creates `Results/<your-name>/`. | Every session. |
| 3. LOAD | Loads `all_graphs.pkl` into `G`, `G_dir`, `B`, `inlets`, `outlets`, `index`. Loads the unweighted graphs by default. For the weighted ones, change `GRAPHS` to `'Results/timeseries_weighted/graphs/all_graphs.pkl'`. | Every session. |
| 4. OVERVIEW | Table of nodes, links, cycles, inlets and outlets per year, to check that loading worked. | Optional. |
| 5. HELPERS | `simple(g)`: simple graph keeping the shortest of parallel links. `pos_km(g)`: positions in km. `draw(g, ax)`: plots on real coordinates. `YEARS`: graph name → year. | Every session. |
| 6. EXAMPLE | Clustering and global efficiency per year, saved as CSV and PNG in your Results folder. | Replace with your work. |
| 7. GIT | Connects to GitHub with your token, pulls, and stages only `team/<your-name>` and `Results/<your-name>`. Commit and push are commented out. | When saving to GitHub (see below). |

Because cell 1 moves into `Mangrove_Maze`, all paths (`Data/`, `Results/`) are relative to the project folder, wherever your notebook is saved.

---

## The shared graphs (`main.ipynb`)

For each image in `Data/` (2016–2026, skipping 2018), `main.ipynb`:
1. resamples every image onto **one common grid** and uses **one common shoreline** (where at least 5 of the 10 water masks agree), so outlets don't move with the tide;
2. extracts the cleaned RivGraph network **G1 … G10**: an undirected node–edge graph, plus a directed copy;
3. builds the baseline **B1 … B10**: a random tree on the same year's nodes (outlets are roots, inlets are leaves);
4. saves everything to `Results/timeseries/graphs/` and plots the parameters over the years against 1,000 baseline trees per year.

### Saved files and formats

Everything is in `Results/timeseries/graphs/`. The same graphs are saved in two formats:

| File | Format | Contents |
|---|---|---|
| `all_graphs.pkl` | **Python pickle** of networkx objects. **Use this in Colab.** | A dict with these keys:<br>• `G`: `{'G1': nx.MultiGraph, …, 'G10': …}`, the real networks, undirected<br>• `G_directed`: the same networks with RivGraph's flow directions (MultiDiGraph; these contain some direction loops)<br>• `B`, `B_directed`: the baseline random trees `{'B1': …, …}`, undirected, and directed child → parent<br>• `inlets`, `outlets`: `{'G1': [node ids], …}`<br>• `index`: a DataFrame with columns `graph`, `baseline`, `date`, `file` |
| `G{k}_{date}.graphml` | **GraphML** (XML). Use this outside Python: Gephi, R/igraph, Cytoscape. | One real network. Nodes: `kind`, `x`, `y`. Links: `length_m`, `width_m`. |
| `B{k}_{date}.graphml` | GraphML | One baseline tree, on the same nodes. Its links may have no length or width. |
| `index.csv` | CSV | Which G/B number goes with which date and image file. |

**Attributes of the real graphs `G`** (in the pickle):

| Where | Attribute | Meaning |
|---|---|---|
| node | `kind` | `'inlet'`, `'outlet'` or `'junction'` |
| node | `pos` | `(x, y)` in metres, UTM zone 45N (EPSG:32645). The GraphML files store these as `x` and `y`. |
| link | `length` | Channel length along the centreline, in metres (`length_m` in GraphML) |
| link | `width` | Mean channel width, in metres (`width_m` in GraphML) |

**Things to know:**
- **The graphs are MultiGraphs.** Two nodes can be joined by more than one channel, and those parallel links are real. Many networkx functions (clustering, some centralities) need a simple graph, so use `simple(g)` from the template. It keeps the shortest parallel link.
- **The pickle needs matching package versions.** It loads reliably only with the same networkx version. Installing from the wheelhouse (template cell 1) guarantees that.
- **GraphML changes node ids.** `nx.read_graphml(path, force_multigraph=True)` gives string node ids (`'12'`, not `12`). Use `node_type=int` to get integers back.
- **These files exist only after a run.** If `Results/timeseries/graphs/` is empty, someone needs to run `main.ipynb` first.

**Load them** (the template's LOAD cell does this):

```python
import pickle
with open('Results/timeseries/graphs/all_graphs.pkl', 'rb') as fh:
    data = pickle.load(fh)
G, B, index = data['G'], data['B'], data['index']
G1 = G['G1']                                               # 2016-10-17

# or one year from GraphML:
import networkx as nx
G1 = nx.read_graphml('Results/timeseries/graphs/G1_2016-10-17.graphml', node_type=int, force_multigraph=True)
```

---

## The weighted graphs (`main_weighted.ipynb`)

This notebook redoes `main.ipynb` on **simplified, weighted** graphs. It doesn't touch `main.ipynb` or `Results/timeseries/`. For each image (same images and the same 2018 skip), it:
1. extracts the network exactly as `main.ipynb` does (same common grid and shoreline);
2. adds weights to each link:
   - `conductance` = width / length, so resistance is L/w: short, wide channels conduct easily;
   - `weight` = conductance ÷ that graph's largest conductance;
3. contracts redundant links, which gives **G1 … G10**:
   - pass 1 contracts the highest-weight outliers, above median + 2 × 1.4826 × MAD, at most 7 per graph;
   - pass 2 contracts links shorter than 120 m;
   - inlets and outlets are never contracted, real loops are kept, and no node has a self-edge;
4. builds baseline trees **B1 … B10**. Each baseline link gets the width and length of a randomly chosen real link from that year;
5. compares both over the years: the unweighted parameters on unweighted copies, and the weighted parameters using `conductance`.

**Use `conductance` for comparisons between years, not `weight`.** `weight` is scaled by each graph's own maximum, so the same value means different things in different years.

### Saved files and formats

Everything is in `Results/timeseries_weighted/graphs/`:

| File | Format | Contents |
|---|---|---|
| `all_graphs.pkl` | Python pickle. **Use this in Colab.** | Same keys as the unweighted pickle: `G`, `G_directed`, `B`, `B_directed`, `inlets`, `outlets`, `index`. `G` and `B` are the simplified, weighted graphs. `G_directed` has the same contractions. `B_directed` has no weights. |
| `G{k}_{date}_weighted.graphml` | GraphML | Simplified graph. Nodes: `kind`, `x`, `y`, `merged`. Links: `weight`, `conductance`, `length_m`, `width_m`. |
| `G{k}_{date}_unweighted.graphml` | GraphML | The same graph, with no link attributes. |
| `B{k}_{date}_weighted.graphml`, `B{k}_{date}_unweighted.graphml` | GraphML | The baseline tree, with and without weights. |
| `index.csv` | CSV | Which G/B number goes with which date and image file. |

`Results/timeseries_weighted/` also holds:
- `graph_counts.csv`: nodes, links and cycles before and after simplification, and how many links were contracted;
- the metric tables, z-scores, trend tests and plots.

**Attributes of the weighted graphs** (in the pickle):

| Where | Attribute | Meaning |
|---|---|---|
| node | `kind`, `pos` | As in the unweighted graphs |
| node | `merged` | `True` if the node was made by contracting links |
| link | `length`, `width` | Metres, as in the unweighted graphs |
| link | `conductance` | width / length. **Use this for weighted measures.** |
| link | `weight` | conductance ÷ that graph's maximum. Only for drawing or ranking within one graph. |

**Load them:**

```python
import pickle
with open('Results/timeseries_weighted/graphs/all_graphs.pkl', 'rb') as fh:
    data_w = pickle.load(fh)
Gw = data_w['G']                                           # {'G1': weighted, simplified MultiGraph, ...}

# or with the helper in sundarbans_simplified_network.py:
from sundarbans_simplified_network import load_graphs, to_unweighted
weighted   = load_graphs('Results/timeseries_weighted/graphs/all_graphs.pkl')                    # with weights
unweighted = load_graphs('Results/timeseries_weighted/graphs/all_graphs.pkl', weighted=False)    # same graphs, no weights
```

For distances in networkx (shortest paths, betweenness, efficiency), use resistance `length / width` as the distance, not `conductance`. For example, set `d['res'] = d['length'] / d['width']` on each link and pass `weight='res'`.

---

## Saving your work: Drive **and** GitHub

**Drive.** Colab saves your notebook to Drive automatically (*File → Save*, or Ctrl+S). Files your code writes to `Results/` land on Drive straight away.

**GitHub.** You must commit and push all your changes to GitHub as well. Do it at the end of every working session, and whenever you finish a piece of analysis.

1. **Use the GIT cell at the end of your notebook** (from the template). It uses your own name, email and `GITHUB_TOKEN`. Don't use the GITHUB cell in `Mangroves.ipynb`: it has Kirthi's details.
2. **Commit only your own files.** Everyone shares one working copy on Drive. The GIT cell stages only `team/<your-name>` and `Results/<your-name>`. Never use `git add -A`, which would also commit teammates' unfinished work.
   - Run the cell once and read the `git status` list.
   - If the list is right, remove the `#` in front of the commit and push lines, write a message, and run it again:

   ```python
   !git commit -m "<your-name>: short description of what changed"
   !git push origin main
   ```

3. **One person at a time.** Everyone shares the same `.git` folder. If git reports `index.lock` or `Another git process seems to be running`, someone else is committing. Wait a minute and try again. Don't delete the lock file while they are working.

### Rules
- Never commit tokens, passwords, or outputs that contain them. Look through `git status` before every commit.
- Don't commit large files that can be rebuilt, such as `wheels_py313/` (already ignored). Ask before adding new images to `Data/`.
- If you change the shared graphs (`main.ipynb`, `main_weighted.ipynb`) or a shared helper (`sundarbans_rg.py`, `tree_nulls.py`, `sundarbans_simplified_network.py`, `mask_tools.py`), say so in the commit message and tell the team. Everyone's results depend on them.

---

## Quick start

1. Do the first-time setup above.
2. Copy `team/TEMPLATE.ipynb` to `team/<your-name>/<your-name>_<topic>.ipynb`, and set `ME` and your GitHub details.
3. Run cells 1–5 to set up and load the shared graphs. For the weighted graphs, point `GRAPHS` at `Results/timeseries_weighted/graphs/all_graphs.pkl`.
4. Replace the EXAMPLE cell with your analysis, saving outputs to `Results/<your-name>/`.
5. Commit and push your folders with the GIT cell at the end of the session.
