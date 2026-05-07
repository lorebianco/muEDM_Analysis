import struct
import numpy as np

DTQ_SPECT = 0x01
DTQ_TIMING = 0x02
DTQ_TSPECT = 0x03
DTQ_COUNT = 0x04
DTQ_SERVICE = 0x2F


class Spectroscopy:
    # Header: uint32, double, double, uint64, uint64,  
    header_fmt = "<IddQQQ"   # I = uint32, d = double, H = uint16
    
    def __init__(self, raw):
        offset = 0    
        # ----- Parse Header -----
        header_size = struct.calcsize(self.header_fmt)
        self.DataQualifier, self.tstamp, self.reltstamp, self.trgID, self.chmask, self.qdmask = struct.unpack(
            self.header_fmt,
            raw[offset : offset + header_size]
        )
        offset += header_size

        # ----- Parse variable-length hit data -----
        charge_list = []
        time_list = []

        for ch in range(64):
            if (self.chmask >> ch) & 1:
                # energy LG/HG
                LG = struct.unpack_from("<H", raw, offset)[0]
                offset += 2
                HG = struct.unpack_from("<H", raw, offset)[0]
                offset += 2
                charge_list.append((LG, HG))

                if self.DataQualifier == 3:
                    ToA = struct.unpack_from("<I", raw, offset)[0]
                    offset += 4
                    ToT = struct.unpack_from("<H", raw, offset)[0]
                    offset += 2
                    time_list.append((ToA, ToT))

        # Convert to numpy structured arrays
        self.charge = np.array(
            charge_list,
            dtype=[("energyLG", np.uint16),
                   ("energyHG", np.uint16)]
        )

        if self.DataQualifier == 3:
            self.time = np.array(
                time_list,
                dtype=[("ToA", np.uint32),
                       ("ToT", np.uint16)]
            )
        else:
            self.time = None

    def __repr__(self):
        return (f"Spectroscopy(DataQualifier={self.DataQualifier}, "
                f"tstamp={self.tstamp}, "
                f"trgID={self.trgID}, "
                f"chmask=0x{self.chmask:016X})")
#end of Spectroscopy


class Timing:
    # Header: uint32, double, uint16
    header_fmt = "<IdH"   # I = uint32, d = double, H = uint16
    
    def __init__(self, raw):
        offset = 0
        
        # ----- Parse Header -----
        header_size = struct.calcsize(self.header_fmt)
        self.DataQualifier, self.fine_tstamp, self.nhits = struct.unpack(
            self.header_fmt,
            raw[offset : offset + header_size]
        )
        offset += header_size

        # ----- Parse variable-length hit data -----
        hits = []
        for _ in range(self.nhits):
            # channel (uint8)
            channel = struct.unpack_from("<B", raw, offset)[0]
            offset += 1

            # tstamp (uint32)
            tstamp = struct.unpack_from("<I", raw, offset)[0]
            offset += 4

            # ToT (uint16)
            tot = struct.unpack_from("<H", raw, offset)[0]
            offset += 2

            hits.append((channel, tstamp, tot))

        # Convert to numpy structured array (optional but useful)
        self.hits = np.array(
            hits,
            dtype=[("channel", np.uint8),
                   ("tstamp", np.uint32),
                   ("tot", np.uint16)]
        )

    def __repr__(self):
        return (f"Timing(DataQualifier={self.DataQualifier}, "
                f"fine_tstamp={self.fine_tstamp}, "
                f"nhits={self.nhits}, hits={self.hits})")
#end of Timing




class Counting:
        # Header: uint32, double, double, uint64, uint32, uint32, uint64  
    header_fmt = "<IddQIIQ"   # I = uint32, d = double, H = uint16
    
    def __init__(self, raw):
        offset = 0    
        # ----- Parse Header -----
        header_size = struct.calcsize(self.header_fmt)
        self.DataQualifier, self.tstamp, self.reltstamp, self.trgID, self.t_or, self.q_or, self.chmask = struct.unpack(
            self.header_fmt,
            raw[offset : offset + header_size]
        )
        offset += header_size

        # ----- Parse variable-length hit data -----
        counts = []
        for i in range(64):
            if (self.chmask>>i) & 1:
                # count per channel (uint32)
                ch_count = struct.unpack_from("<I", raw, offset)[0]
                offset += 4
                counts.append(ch_count)

        # Convert to numpy structured array (optional but useful)
        self.counts = np.array(counts, dtype=np.uint32)
        

    def __repr__(self):
        return (f"Counting(DataQualifier={self.DataQualifier}, "
                f"tstamp={self.tstamp}, "
                f"trgID={self.trgID}, "
                f"chmask=0x{self.chmask:016X})")
#end of Counting


class Service:
    # d = double, Q = uint64, H = uint16, B = uint8, I = uint32, f = float
    # ch_trg_cnt = array di FERSLIB_MAX_NCH_5202 uint32
    header_fmt = f"<dQHBB64I" \
                 "III" + "f" * 8 + "BBBB" + "HH" + "2Q" + "III" * 3

    def __init__(self, raw):
        offset = 0
        # ----- Parse fixed fields -----
        # tstamp_us, update_time, pkt_size, version, format
        self.tstamp_us, self.update_time, self.pkt_size, self.version, self.format = struct.unpack_from("<dQHBB", raw, offset)
        offset += struct.calcsize("<dQHBB")

        # ch_trg_cnt array
        self.ch_trg_cnt = np.array(
            struct.unpack_from(f"<64I", raw, offset),
            dtype=np.uint32
        )
        offset += struct.calcsize(f"<64I")

        # q_or_cnt, t_or_cnt, TotTrg_cnt
        self.q_or_cnt, self.t_or_cnt, self.TotTrg_cnt = struct.unpack_from("<III", raw, offset)
        offset += struct.calcsize("<III")

        # temperature fields: tempFPGA, tempBoard, tempTDC[2], tempHV, tempDetector, hv_Vmon, hv_Imon
        self.tempFPGA, self.tempBoard, tdc0, tdc1, self.tempHV, self.tempDetector, self.hv_Vmon, self.hv_Imon = struct.unpack_from("<ff2fff f f", raw, offset)
        self.tempTDC = np.array([tdc0, tdc1], dtype=np.float32)
        offset += struct.calcsize("<ff2fff f f")

        # hv_status fields: hv_status_on, hv_status_ramp, hv_status_ovv, hv_status_ovc
        self.hv_status_on, self.hv_status_ramp, self.hv_status_ovv, self.hv_status_ovc = struct.unpack_from("<BBBB", raw, offset)
        offset += struct.calcsize("<BBBB")

        # Status, TDCROStatus
        self.Status, self.TDCROStatus = struct.unpack_from("<HH", raw, offset)
        offset += struct.calcsize("<HH")

        # ChAlmFullFlags[2]
        self.ChAlmFullFlags = np.array(struct.unpack_from("<2Q", raw, offset), dtype=np.uint64)
        offset += struct.calcsize("<2Q")

        # ReadoutFlags, TotTrg_cnt, RejTrg_cnt, SupprTrg_cnt
        self.ReadoutFlags, self.TotTrg_cnt, self.RejTrg_cnt, self.SupprTrg_cnt = struct.unpack_from("<IIII", raw, offset)
        offset += struct.calcsize("<IIII")

    def __repr__(self):
        return (f"Service(tstamp_us={self.tstamp_us}, "
                f"update_time={self.update_time}, "
                f"pkt_size={self.pkt_size}, "
                f"version={self.version}, "
                f"format={self.format})")
#end of Service
