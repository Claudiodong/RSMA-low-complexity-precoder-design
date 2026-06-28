# RSMA Repository Reorganization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reorganize the RSMA MATLAB project into a clearer academic repository and replace the minimal README with an illustrated project guide.

**Architecture:** Keep each original figure experiment self-contained so MATLAB scripts can still run from their experiment folder. Move rendered images into `figures/`, move MATLAB experiment folders into `experiments/`, and keep supporting planning documentation under `docs/`.

**Tech Stack:** MATLAB, MAT result files, GitHub Markdown, existing PNG outputs.

---

### Task 1: Reorganize Project Folders

**Files:**
- Move: `figure 4.1/` to `experiments/figure-4.1-snr-initialisation/`
- Move: `figure 4.2/` to `experiments/figure-4.2-ao-convergence/`
- Move: `figure 4.3 4.4 4.6/` to `experiments/figures-4.3-4.4-4.6-low-complexity/`
- Move: `figure 4.7/` to `experiments/figure-4.7-adaptive-common-streams/`
- Move: `figure 4.8 4.9/` to `experiments/figures-4.8-4.9-power-distribution/`
- Move: `figure 4.10/` to `experiments/figure-4.10-large-scale-comparison/`
- Move: `png/` to `figures/`

- [x] **Step 1: Create destination directories**

Run: `mkdir -p experiments figures`

Expected: `experiments/` and `figures/` exist.

- [x] **Step 2: Move the original folders with git history**

Run each `git mv` command listed in the Files section.

Expected: `git status --short` shows renames rather than deleting the project content.

### Task 2: Replace README

**Files:**
- Modify: `README.md`

- [x] **Step 1: Replace the minimal README**

Write a GitHub Markdown README with:
- project title and MSc context
- short explanation of RSMA and the low-complexity contribution
- embedded figures from `figures/`
- repository map for each experiment folder
- reproduction instructions for MATLAB and CVX
- notes on results and academic references

- [x] **Step 2: Check README paths**

Run: `rg "experiments/|figures/" README.md`

Expected: all linked project paths point to existing files or folders.

### Task 3: Verify and Publish

**Files:**
- Inspect: all renamed folders and `README.md`

- [x] **Step 1: Verify repository status**

Run: `git status --short`

Expected: folder renames, README update, and this plan file are visible.

- [x] **Step 2: Commit changes**

Run: `git add -A && git commit -m "docs: reorganize RSMA project presentation"`

Expected: commit succeeds.

- [x] **Step 3: Push to GitHub**

Run: `git push origin main`

Expected: GitHub main branch receives the reorganization commit.
