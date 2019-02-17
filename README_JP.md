# elfin v2 [![Build Status](https://travis-ci.com/joy13975/elfin.svg?branch=master)](https://travis-ci.com/joy13975/elfin)

Elfinはリピートタンパク質の組み合わせにより、より大きいタンパク質をデザインすることに補助するツールセットです。 OS、ベンダやアーチ固有の仮定をしていませんが、Linux、MacOS、およびWindows no WSL上でのみテストされています。

Elfinは４つのリポジトリに割れています：
 1. [elfin-ui](https://github.com/joy13975/elfin-ui) - GUI。
 2. [elfin-library](https://github.com/joy13975/elfin-library) - 公開なデータベース。
 3. [elfin-data](https://github.com/joy13975/elfin-data) - 非公開なデータベース。
 4. [elfin-solver](https://github.com/joy13975/elfin-solver) - 自動デザインの機械学習アルゴリズム。

このメインのリポジトリはデータ処理のスクリプトを含んでいます。


# elfinって何？

Elfinはタンパク質を築くための補助ツールセットです。Elfinの特色は、機械学習の自動デザインアルゴリズムが実装されているというところです。Elfinの中に使われている論理は[この学術論文](https://www.sciencedirect.com/science/article/pii/S1047847717301417)に記載しております。

[スキップ: 使い方法](#使い方法)

要に、リピートタンパク質は他のリピートタンパク質と接続して1つの大きいタンパク質を築くことができます。Elfinはこの理論を利用し、リピートの組み合わせによりできるだけ使用者がGUIで描いた形と近い形を作成します。

下記の方々のおかげでこのプロジェクトが成し遂げられました：
 - スーパーバイザー：Fabio Parmeggiani (UoB)
 - タンパク質データ研究者＆提供者：TJ Brunette (UoW), David Baker (UoW)
 - テクニカルコンサルタント：Simon McIntosh-Smith (UoB)

(UoB: University of Bristol, UoW: University of Washington)

Elfin v2は概念証明型のv1([branch v1](https://github.com/joy13975/elfin/tree/v1))より複雑なデザインができる進化型です。

データベースにある新しいリピートタンパク質の学術論文は審査中のため、データの提供者の依頼でタンパク質のPDBファイルは[非公開のリポジトリ](https://github.com/joy13975/elfin-db)にアップロードしております。

必要であれば、[Fabio Parmeggiani](https://github.com/parmef)さん（fabio.parmeggiani@bristol.ac.uk）にメールをお送りください。PDBファイルがなくてもElfinは動けるのです。但し、最後のステージにタンパク質のPDB或いはCIFファイルを出力することはできません。

![alt tag](resources/diagrams/ProteinBristol.png)
Figure 1: 手で書いた「Bristol」の形をタージェットに、elfin (v1)で自動的にデザインさせたタンパク質。可視化は[PyMol](https://pymol.org)で実現しました.

## プロジェクト現状

Elfinバージョン２の機能は本とんど完全ですが、いくつかのノンクリティカルTODOはそれぞれのリポジトリのIssuesに記録しております。

## 使い方法
### 目次
   1. [前提ソフトウェア](#1-前提ソフトウェア)
   2. [セットアップ作業](#2-セットアップ作業)
   3. [デザインのワークフロー](#3-デザインのワークフロー)
   4. [スクリプティング](#4-スクリプティング)

## 1. 前提ソフトウェア
まず、下記のソフトウェアをインストールしてください：
#### 必須
1. [Python 3+](https://www.python.org/downloads/)
2. [Virtualenv](https://virtualenv.pypa.io/en/stable/)
3. [Blender](https://www.blender.org/)
4. [gcc-5+](https://gcc.gnu.org/)

#### 任意
1. [PyMOL](https://www.pymol.org) - タンパク質の可視化のため。
2. [Rosetta](https://www.rosettacommons.org/software/license-and-download) - タンパク質のデータの最適化のため.

## 2. セットアップ作業

下記のコマンドを実行してください:
```Bash
bash <(curl -s https://raw.githubusercontent.com/joy13975/elfin/master/setup_elfin)
```

ご注意：このスクリプトは自動的に非公開の[elfin-data](https://github.com/joy13975/elfin-data)からデータをダウンロードするため、Githubのユーザー名とパスワードをご入力をいただきます。もしelfin-dataのアクセスの許可がなければ、スキップ（Enterを二回）してもかまわないです。

## 3. デザインのワークフロー

#### 立体幾何形を描く
始めに、使用者がタンパク質で築きたい形を[elfin-ui](https://github.com/joy13975/elfin-ui)で描きます。ドキュメンテーションはそのリポジトリに含まれています。

その後、JSONファイルを出力します（elfin-uiコマンド: #exp）。

#### 自動的にデザイン
そして、[elfin-solver](https://github.com/joy13975/elfin-solver)を使い、求めている立体幾何形を自動的にデザインさせます。elfin-solverのデザイン結果はもう一つのJSONファイルに出力します。

#### デザインを直す
Blenderの中で、新しいファイルを開き、 先ほどelfin-solverの結果JSONファイルを入力します (elfin-uiコマンド: #imp)。

必要であれば、デザインを直します(おそらく隙をなくすやネットワークを接ぎます)。

直したデザインをもう一度出力する(#exp)上、今回はPath Guideオブジェクトのないようにしてください。 次のステージはPath Guideオブジェクトを含んでいるJSONファイルを受け取りません。

#### PDB/CIFの出力
最後に、下記のコマンドを実行してください：

```
. ./activate  # Activates venv
stitch.py <PATH_TO_YOUR_EXPORTED_SOLUTION_JSON>
```

ご注意：```. ./activate```は`venv`に切り替えていない場合のみ必要です。

#### データ前処理

タンパク質のデータは既に前処理しており、[elfin-data](https://github.com/joy13975/elfin-data)にアップロードしておりますので、ほとんどの方にとって前処理は必要がありません。とはいえ、新しいタンパク質が加えられたり、前処理の仕方が変わったりすれば、前処理をやり直すは必要となります。 

そういう場合には[elfin-data](https://github.com/joy13975/elfin-data)(private)へお越しください。

### 4. スクリプティング

Elfinのデータ処理クラスはあなたのスクリプトやPythonのインタラクティブシェルでの利用が可能です。 

改めて、virtualenvの環境を確かめましょう:
```
. ./activate
```

下記はStitcherを使用する例です：

```Python
import elfinpy.utilities as utils
import elfinpy.pdb_utilities as pdb_utils
import elfinpy.stitch as stitch

struct = stitch.Stitcher(
    spec=utils.read_json('resources/examples/half_snake_2x1h_deposit_test.json'),
    xdb=utils.read_json('resources/xdb.json'),
    pdb_dir='./resources/pdb_aligned/',
    cappings_dir='./resources/pdb_relaxed/cappings',
    metadata_dir='./resources/metadata/',
    show_fusion=False,
    disable_capping=False,
    skip_unused=False).run()

as_cif = True

if as_cif:
    pdb_utils.save_cif(struct=struct, path='test.cif')
else:
    cid = 'A'
    for c in struct[0].get_chains():
        if cid > 'Z':
            raise ValueError('Too many chains for PDB format')
            # You might want to extend number of chains by using lowercase
            # alphabets.

        c.id = cid
        cid = chr(ord(cid) + 1)

    pdb_utils.save_pdb(struct=struct, path='test.pdb')
```