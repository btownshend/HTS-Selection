from rdkit.ML.InfoTheory import InfoBitRanker


def getFragRanks(fcat, fps, acts, ntop=5):
    ranker = InfoBitRanker(len(fps[0]), 2)
    print("Num of targets active: %.0f, inactive: %.0f" % (sum(acts), len(acts) - sum(acts)))
    for i, fp in enumerate(fps):
        ranker.AccumulateVotes(fp, int(acts[i]))
    top5 = ranker.GetTopN(ntop)
    for fid, gain, n0, n1 in top5:
        ifid = int(fid)
        hitstr = "".join(["+" if fps[i].GetBit(ifid) else "-" for i in range(len(acts)) if acts[i]])
        print('%5d' % ifid, '%.3f' % gain, '%2d' % int(n0), '%2d' % int(n1), hitstr, fcat.GetEntryDescription(int(fid)))

    return top5[0][0]
