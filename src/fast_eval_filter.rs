impl Iterator for RecordsIterator {
    type Item = Result<VariantRecord, bcf::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = self.reader.empty_record();
        match self.reader.read(&mut record) {
            Some(Ok(_)) => {
                let mut info = HashMap::new();
                for tag in record.header().info_tags() {
                    if let Ok(value) = record.info(tag.as_bytes()).string() {
                        if let Some(value) = value {
                            // Assuming the first value for simplicity; needs adjustment for multiple values
                            info.insert(tag, String::from_utf8_lossy(value[0]).to_string());
                        }
                    }
                }
                // FastEval expression filtering integration placeholder
                Some(Ok(VariantRecord { info }))
            },
            Some(Err(e)) => Some(Err(e)),
            None => None,
        }
    }
}- src/fast_eval_filter.rs:
