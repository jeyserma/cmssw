import FWCore.ParameterSet.Config as cms

muonvtxagniso04 = cms.EDProducer("MuonVtxAgnosticIsoProducer",
    muonInputTag = cms.InputTag("linkedObjects","muons"),
    pfCandidateInputTag = cms.InputTag("packedPFCandidates"),
    maxdr = cms.double(0.4),
    mindr = cms.double(0.005),
    maxdeltaz = cms.double(0.2),
    minptchg = cms.double(0.0),
    minptneu = cms.double(0.0),
    minptpho = cms.double(0.0),
    minptpu = cms.double(0.0),
    calculateNeutralPhoton = cms.int32(1)
)

muonvtxagniso03 = muonvtxagniso04.clone()
muonvtxagniso03.maxdr = cms.double(0.3)

muoncrosscheckiso04 = cms.EDProducer("MuonCrossCheckIsoProducer",
    muonInputTag = muonvtxagniso04.muonInputTag,
    pfCandidateInputTag = muonvtxagniso04.pfCandidateInputTag,
    pvInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
    maxdr = muonvtxagniso04.maxdr,
    mindr = muonvtxagniso04.mindr,
    maxdeltaz = muonvtxagniso04.maxdeltaz,
    minptchg = muonvtxagniso04.minptchg,
    minptneu = muonvtxagniso04.minptneu,
    minptpho = muonvtxagniso04.minptpho,
    minptpu = muonvtxagniso04.minptpu,
    calculateNeutralPhoton = muonvtxagniso04.calculateNeutralPhoton
)

muoncrosscheckiso03 = muoncrosscheckiso04.clone()
muoncrosscheckiso03.maxdr = muonvtxagniso03.maxdr
