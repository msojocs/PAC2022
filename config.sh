#!/bin/sh
# using sh config.sh
# default commit id
commitID=92d25e35d9bf1a6b16f7d0758f25d48ace11e5b9
if [ "$1" != "" ]; then
    commitID="$1"
fi
echo "commit id: ${commitID}"
if [ ! -f "/tmp/vscode-server/${commitID}/vscode-server-linux-x64.tar.gz" ]; then
    echo "can not find extension file !"
    exit 1
fi
cp "/tmp/vscode-server/${commitID}/vscode-server-linux-x64.tar.gz" ~/.vscode-server/bin/vscode-server.tar.gz
cd ~/.vscode-server/bin
tar -zxf vscode-server.tar.gz
rm -rf "${commitID}"
mv vscode-server-linux-x64 "${commitID}"

rm -rf ~/.vscode-server/extensions
cp "/tmp/vscode-server/extensions/extensions.tar.gz" ~/.vscode-server/extensions.tar.gz
cd ~/.vscode-server
tar -izxf extensions.tar.gz
ls -l ~/.vscode-server/extensions

echo "done"
ls -l ~/.vscode-server/bin