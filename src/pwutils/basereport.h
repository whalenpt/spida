// basereport.h
#pragma once

#include <cassert>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <string_view>
#include "pwutils/pwdefs.h"
#include "pwutils/reporthelper.h"

namespace pw{

enum class TrackType { Max, Min };

class ReportBase{
    public:
        explicit ReportBase(const std::string& nm)
          : m_name(nm) {}
        virtual ~ReportBase() = default;

        void report(std::ofstream& os) const;
        void report(std::ofstream& os, unsigned rep_num) const;
        friend std::ofstream& operator<<(std::ofstream& os, const ReportBase& def);
        friend std::ofstream& operator<<(std::ofstream& os, const ReportBase* def);
        friend std::ofstream& operator<<(std::ofstream& os, const std::unique_ptr<ReportBase>& def);

        void setName(std::string_view nm) {m_name = nm;}
        void setItem(const std::string& key, double val);
        void setItem(const std::string& key, float val);
        void setItem(const std::string& key, int val);
        void setItem(const std::string& key, const std::string&);
        void removeItem(const std::string&);
        void setFileExtension(std::string_view ext) {m_extension = ext;}
        void setDirPath(const std::filesystem::path& dirpath) {
            pw::createDirectory(dirpath, false);
            m_dirpath = dirpath;
        }

        std::string getName() const {return m_name;}
        const metadataMap& getMetadata() const {return m_metadata_map;}
        std::filesystem::path dirpath() const {return m_dirpath;}
        std::filesystem::path path() const {
            return pw::filePath(m_dirpath, m_name, m_extension);
        }
        std::filesystem::path path(int rep_num) const {
            return pw::filePath(m_dirpath, m_name, rep_num, m_extension);
        }

    private:
        std::string m_name;
        metadataMap m_metadata_map;
        std::string m_extension{"json"};
        virtual void reportImplement(std::ofstream& os) const = 0;
        std::filesystem::path m_dirpath{"outfolder"};
};

class ReportData1D : public ReportBase {
    using ReportBase::ReportBase;
};

class ReportData2D : public ReportBase {
    public:
        using ReportBase::ReportBase;
        ~ReportData2D() override = default;
        unsigned getStrideX() const {return m_strideX;}
        unsigned getStrideY() const {return m_strideY;}
        void setStrideX(unsigned s) { assert(s >= 1); m_strideX = s; }
        void setStrideY(unsigned s) { assert(s >= 1); m_strideY = s; }
    private:
        unsigned m_strideX{1};
        unsigned m_strideY{1};
};

class TrackData : public ReportBase {
    public:
        TrackData(const std::string& nm, TrackType ttype)
          : ReportBase(nm), m_ttype(ttype) {}
        ~TrackData() override = default;
        virtual void updateTracker(double t) = 0;
        TrackType getTrackType() const {return m_ttype;}
        void setTrackType(TrackType t) {m_ttype = t;}
    private:
        TrackType m_ttype;
};

}
