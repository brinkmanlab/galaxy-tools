#!/usr/bin/env bash
ln -sf output/`ls -1 output | tail -n1` final
mv "`ls -1 final/*.backbone | tail -n1`" final.backbone
mv "`ls -1 final/*_contigs.tab | tail -n1`" final.contigs.tab
mv "`ls -1 final/*.fas | tail -n1`" final.reordered

FEATURES="`ls -1 final/*_features.tab | tail -n1`"
if [[ -w ${FEATURES} ]]; then
mv  final.features.tab
fi

