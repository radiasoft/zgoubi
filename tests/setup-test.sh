#!/usr/bin/env bash
cd ${zgoubi_root}
git apply tests/prototypical-patch.diff

pushd exemples/KEYWORDS/FIT/FIT_embedded_in_REBELOTE/
  ln -s centeredHelix_FIT_save_nofinal_150226.res zgoubi.dat
popd
