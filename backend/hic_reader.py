import numpy as np

try:
    import hicstraw
    HIC_AVAILABLE = True
except ImportError:
    try:
        import straw as hicstraw
        HIC_AVAILABLE = True
    except ImportError:
        HIC_AVAILABLE = False


class HicFileReader:
    def __init__(self, hic_path):
        if not HIC_AVAILABLE:
            raise ImportError("hic-straw or hicstraw package required for .hic support. Install with: pip install hic-straw")
        self.path = hic_path
        self.hic = hicstraw.HiCFile(hic_path)
        self.resolutions = sorted(self.hic.getResolutions(), reverse=True)
        self._chrom_info = self._load_chroms()

    def _load_chroms(self):
        chroms = self.hic.getChromosomes()
        info = {}
        for c in chroms:
            name = c.name.replace('chr', '').replace('Chr', '')
            if name == 'All' or name == 'MT' or name == 'M':
                continue
            info[name] = {"length": c.length, "name": name}
        return info

    def get_all_chroms(self):
        chrom_info = []
        offset = 0
        for name, ci in self._chrom_info.items():
            chrom_info.append({"name": name, "length": ci["length"], "offset": offset})
            offset += ci["length"]
        return chrom_info

    def get_genome_index(self):
        chrom_list = self.get_all_chroms()
        return {
            "total_length": sum(c["length"] for c in chrom_list),
            "chromosomes": chrom_list
        }

    async def fetch_range_data(self, res, chrom, start, end, start2=None, end2=None):
        if start2 is None:
            start2 = start
        if end2 is None:
            end2 = end
        chrom_name = f"chr{chrom}" if not chrom.startswith("chr") else chrom
        try:
            mzd = self.hic.getMatrixZoomData(chrom_name, chrom_name, "observed", "KR", "BP", res)
            records = mzd.getRecords(start, end, start2, end2)
            results = []
            for rec in records:
                if rec.counts > 0:
                    results.append({"x": rec.binX * res, "y": rec.binY * res, "v": float(rec.counts)})
            return results, res
        except Exception:
            return [], res

    def get_contact_matrix(self, chrom, start, end, resolution):
        chrom_name = f"chr{chrom}" if not chrom.startswith("chr") else chrom
        try:
            mzd = self.hic.getMatrixZoomData(chrom_name, chrom_name, "observed", "KR", "BP", resolution)
            records = mzd.getRecords(start, end, start, end)
            n_bins = (end - start) // resolution
            matrix = np.zeros((n_bins, n_bins))
            for rec in records:
                i = (rec.binX * resolution - start) // resolution
                j = (rec.binY * resolution - start) // resolution
                if 0 <= i < n_bins and 0 <= j < n_bins:
                    matrix[i, j] = rec.counts
                    matrix[j, i] = rec.counts
            return matrix
        except Exception:
            n_bins = (end - start) // resolution
            return np.zeros((n_bins, n_bins))

    def close(self):
        pass
