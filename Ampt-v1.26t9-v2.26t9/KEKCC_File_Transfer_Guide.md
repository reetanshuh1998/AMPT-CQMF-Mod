# KEKCC File Transfer Guide: The Termbin Pipeline

This guide documents the "outbound pull" method for transferring files into secure high-performance computing clusters (like KEKCC) that block inbound `scp`/SSH connections to compute nodes.

## The Concept
Clusters block **inbound** connections (pushing from your laptop to the cluster) but usually allow **outbound** connections (the cluster reaching out to the internet). By uploading your files to a temporary text host (`termbin.com`) and asking the cluster to pull them down via `curl`, you safely bypass the inbound firewall.

---

## Scenario 1: Uploading a Single Text File
Use this for simple files like a Python script, a `.sh` file, or a config file.

**Step 1: On your Local Machine**
Pipe the file's contents to termbin using `nc` (netcat):
```bash
cat my_script.py | nc termbin.com 9999
```
*Output: `https://termbin.com/abcd`*

**Step 2: On KEKCC**
Use `curl` to pull the raw text and save it:
```bash
curl -s https://termbin.com/abcd > my_script.py
```

---

## Scenario 2: Uploading Directories or Binary Files
Use this for folders, `.csv` data, `.pdf`s, compiled binaries, or `.tar.gz` archives. You **must** Base64 encode these files; otherwise, the binary data will be corrupted when passed through the text-only termbin server.

**Step 1: On your Local Machine**
Compress your files into a tarball, encode it to Base64, and pipe it to termbin all in one line:
```bash
tar -czf my_bundle.tar.gz file1.csv file2.txt my_folder/
base64 my_bundle.tar.gz | nc termbin.com 9999
```
*Output: `https://termbin.com/wxyz`*

**Step 2: On KEKCC**
Download the Base64 text, decode it back to binary, and unpack the tarball directly from memory:
```bash
curl -s https://termbin.com/wxyz | base64 -d | tar -xzf -
```
*(The `-` at the end tells tar to extract from the incoming data stream instead of a saved file.)*

---

## Security Warning
> **WARNING:**
> Files uploaded to `termbin.com` are technically public if someone guesses your random 4-letter URL. They automatically expire and are deleted after 30 days. 
> 
> **Do not use this method for:**
> - Private SSH keys
> - Passwords or credentials
> - Sensitive personal data
> 
> It is perfectly safe for research scripts, public data tables, and AMPT code.
