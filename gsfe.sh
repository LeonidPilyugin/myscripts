#!/bin/zsh

source $HOME/.zshrc

mkdir -p $HOME/results/gsfe

file=/tmp/$(hexdump -vn16 -e'4/4 "%08x" 1 "\n"' /dev/urandom).toml
cat <<EOF > "$file"
[task]
id = "gsfe.$1"
comment = "Get GSFE in plane"

[task.backend]
type = "screen"

[task.wrapper]
type = "default"
command = "gsfe.py $HOME/descriptors/gsfe/$1.toml $HOME/results/gsfe/$1.csv"
mixed_stdout = false
EOF

conda activate spmi
python3 "$HOME/opt/spmi/src/spmi/app.py" load "$file"
python3 "$HOME/opt/spmi/src/spmi/app.py" start "gsfe.$1"

conda deactivate

rm "$file"
