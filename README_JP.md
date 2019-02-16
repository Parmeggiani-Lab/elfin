# elfin v2 [![Build Status](https://travis-ci.com/joy13975/elfin.svg?branch=master)](https://travis-ci.com/joy13975/elfin)

Elfinはリピートタンパク質の組み合わせにより、より大きいタンパク質をデザインすることに補助するツールセット。 ElfinはOS、ベンダ、
アーチ固有の仮定をしていませんが、Linux、MacOS、およびWSL上でのみテストされています。

Elfinは４つのリポジトリに割れています：
 1. [elfin-ui](https://github.com/joy13975/elfin-ui) - GUI。
 2. [elfin-library](https://github.com/joy13975/elfin-library) - 公開なデータベース。
 3. [elfin-data](https://github.com/joy13975/elfin-data) - 非公開なデータベース。
 4. [elfin-solver](https://github.com/joy13975/elfin-solver) - 自動デザインの機械学習アルゴリズム。

このメインのリポジトリはデータ処理のスクリプトを含んでいます。

# elfinって何？

Elfinは私の卒業研究プロジェクトから生まれたタンパク質デザインツールセット。Elfinの中に使われている論理は[この学術論文](https://www.sciencedirect.com/science/article/pii/S1047847717301417)に記載しております。

[スキップ: 始めましょう](#2-前提ソフトウェア)

要に、リピートタンパク質は他のリピートタンパク質と接続して1つの大きいタンパク質になることができます。Elfinはこの理論を利用し、リピートの組み合わせによりできるだけ使用者がGUIで描いた形と近い形を作成します。

下記の方々のおかげでこのプロジェクトが成し遂げられました。
 - スーパーバイザー：Fabio Parmeggiani (UoB)
 - タンパク質データ研究者＆提供者：TJ Brunette (UoW), David Baker (UoW)
 - テクニカルコンサルタント：Simon McIntosh-Smith (UoB)

(UoB: University of Bristol, UoW: University of Washington)

Elfin v2は概念証明型のv1([branch v1](https://github.com/joy13975/elfin/tree/v1))より複雑なデザインができる進化型です。

データベースにある新しいリピートタンパク質の学術論文は審査中のため、データの提供者の依頼でタンパク質のPDBファイルは[非公開のリポジトリ](https://github.com/joy13975/elfin-db)にアップロードしております。

必要になれば、[Fabio Parmeggiani](https://github.com/parmef)さんにメールをお送りください。PDBファイルがなくてもElfinは動けるのです。但し、最後のステージにタンパク質のPDB或いはCIFファイルを出力することはできません。

![alt tag](resources/diagrams/ProteinBristol.png)
Figure 1: 手で書いた「Bristol」の形をインプットとして、Elfinが自動的にデザインしたタンパク質。 可視化は[PyMol](https://pymol.org)で実現しました.

### 目次
1. [プロジェクト現状](#1-プロジェクト現状)

2. [前提ソフトウェア](#2-前提ソフトウェア)

3. [インストール](#3-インストール)

4. [タンパク質デザインGUI](#4-タンパク質デザインGUI)

5. [機械学習の自動的デザインアルゴリズム](#5-機械学習の自動的デザインアルゴリズム)

6. [タンパク質の全処理](#6-タンパク質の全処理)

7. [デザイン結果の出力](#7-デザイン結果の出力)

## 1. プロジェクト現状

Elfinバージョン２の機能は本とんど完全ですが、いくつかのノンクリティカルTODOはそれぞれのリポジトリのIssuesに記録しております。

## 2. 前提ソフトウェア
#### 必須
1. [Python 3+](https://www.python.org/downloads/)
2. [Virtualenv](https://virtualenv.pypa.io/en/stable/)
3. [Blender](https://www.blender.org/)
4. [gcc-5+](https://gcc.gnu.org/)

#### 任意
1. [PyMOL](https://www.pymol.org) - タンパク質の可視化のため。
2. [Rosetta](https://www.rosettacommons.org/software/license-and-download) - タンパク質のデータの最適化のため.

## 3. インストール

下記のコマンドを実行してください:
```Bash
bash <(curl -s https://raw.githubusercontent.com/joy13975/elfin/master/setup_elfin)
```

ご注意：このスクリプトは自動的に非公開の[elfin-data](https://github.com/joy13975/elfin-data)からデータをダウンロードするため、Githubのユーザー名とパスワードをご入力をいただきます。もしelfin-dataのアクセスの許可がなければ、スキップ（Enterを二回）してもかまわないです。

## 4. タンパク質デザインGUI

[elfin-ui](https://github.com/joy13975/elfin-ui)へお越しください。

## 5. 機械学習の自動的デザインアルゴリズム

[elfin-solver](https://github.com/joy13975/elfin-solver)へお越しください。

## 6. タンパク質の全処理

（非公開）[elfin-data](https://github.com/joy13975/elfin-data)へお越しください。

## 7. デザイン結果の出力

If you want to invoke `stitch.py` on a v1 elfin output json file, first convert it to a `stitch.py` readable format. Taking `resources/examples/horns_output.json` as an example. Run at elfin root:
```Bash
. ./activate
./elfinpy/v1_design_convert.py resources/examples/horns_output.json
```

And then you will be able to stitch the v2 output:
```Bash
./elfinpy/stitch.py resources/examples/horns_output.v2.json
```

If you want to process solver output entirely in code, below is an example:


```Bash
. ./activate # make sure you're in elfin root and in virtual env
```

Then either in a script or an interactive shell:

```Python
import elfinpy.v1_design_convert as converter
import elfinpy.utilities as utils
import elfinpy.pdb_utilities as pdb_utils
import elfinpy.stitch as stitch

# Normally, you'd need to read a JSON file for the output data from solver.
# Alternatively, create the solver output in Python.

input_json = utils.read_json('resources/examples/horns_output.json')
graphs = converter.v1_to_v2(input_json, 'resources/xdb.json')

# Convert List[graph] to dictionary
graphs_dict = utils.to_dict(graphs)

# Run stitch synth
syn = stitch.Synthesiser(graphs_dict, 
    pdb_dir='./resources/pdb_aligned/',
    cappings_dir='./resources/pdb_relaxed/cappings',
    metadata_dir='./resources/metadata/',
    show_fusion=False,
    disable_capping=False)
struct = syn.run()

as_cif = True

if as_cif:
    pdb_utils.save_cif(struct=struct, path='test.cif')
else:
    # To save as PDB, the chain ID needs to be changed to
    # one that PDB format accepts i.e. a single character.
    # The following line assumes there's only 1 chain. By
    # right in a v1 format there should only be 1 chain 
    # anyway.
    [c for c in struct[0].get_chains()][0].id = 'a'
    pdb_utils.save_pdb(struct=struct, path='test.pdb')
```