"""Exercise matched-site QC against real tabix files, including boundaries."""
import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
import pysam

SOURCE=Path(__file__).resolve().parents[1]/'workflow/scripts/quality_control/conversion_qc.py'
spec=importlib.util.spec_from_file_location('conversion_qc',SOURCE)
qc=importlib.util.module_from_spec(spec);spec.loader.exec_module(qc)


class ConversionQC(unittest.TestCase):
    def test_shared_sites_units_pc_and_reproducibility(self):
        with tempfile.TemporaryDirectory(prefix='foodie_qc_test_') as tmp:
            root=Path(tmp)
            bed=root/'open.bed';bed.write_text('chr1\t100\t110\nchr1\t105\t112\n')
            paths=[]
            for sample in range(3):
                raw=root/f's{sample}.txt'
                lines=[]
                for pos in range(99,114):
                    depth=1 if sample==2 and pos==105 else 20
                    converted=(sample+1)*4 if depth==20 else 1
                    lines.append(f'chr1\t{pos}\t{pos+1}\t{converted}\t{depth-converted}\t{converted/depth}\tC\n')
                raw.write_text(''.join(lines))
                gz=str(raw)+'.gz';pysam.tabix_compress(str(raw),gz)
                # Match existing workflow's point-index coordinate convention.
                pysam.tabix_index(gz,seq_col=0,start_col=1,end_col=1,zerobased=False)
                paths.append(gz)
            info=dict(pc_control=paths[0],open_regions=str(bed),samples=[dict(sample_id='dose1',concentration_M=1e-6,path=paths[1]),dict(sample_id='dose2',concentration_M=2e-6,path=paths[2])])
            specification=dict(options=dict(n_sites=100,min_depth=10,seed=42,sample_regions=100),tfs={'TF':info})
            first=qc.run_qc(specification,root/'first')
            second=qc.run_qc(specification,root/'second')
            sites=pd.read_csv(root/'first/TF.sampled_sites.tsv',sep='\t')
            self.assertEqual(set(sites.start),set(range(100,112))-{105})
            self.assertEqual(set(first.n_sites),{11})
            self.assertTrue(np.allclose(first.concentration_uM,[0,1,2]))
            self.assertTrue(np.allclose(first.mean_conversion,[.2,.4,.6]))
            self.assertEqual(len(sites),len(sites[['chrom','start']].drop_duplicates()))
            pd.testing.assert_frame_equal(first,second)
            pd.testing.assert_frame_equal(sites,pd.read_csv(root/'second/TF.sampled_sites.tsv',sep='\t'))
            invalid=copy.deepcopy(specification);invalid['options']['min_depth']=1000
            with self.assertRaisesRegex(ValueError,'No shared'):
                qc.run_qc(invalid,root/'empty')
            self.assertFalse((root/'empty/sampling_metadata.json').exists())


if __name__=='__main__':
    unittest.main()
