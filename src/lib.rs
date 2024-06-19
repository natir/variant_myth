//! A variant annotator.

#![warn(missing_docs)]

/* std use */

/* crate use */

/* project use */

/* mod declaration */
pub mod annotation;
pub mod annotations_db;
pub mod cli;
pub mod error;
pub mod interval_tree;
pub mod variant;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
