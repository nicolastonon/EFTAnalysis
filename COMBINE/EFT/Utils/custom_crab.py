def custom_crab(config):
    print '>> Customizing the crab config'
    from CRABClient.UserUtilities import getUsernameFromSiteDB
    config.Site.storageSite = 'T2_DE_DESY'
    config.Data.outLFNDirBase = '/store/user/%s/EFT/' % (getUsernameFromSiteDB())

